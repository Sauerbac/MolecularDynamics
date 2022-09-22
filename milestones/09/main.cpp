#include "Atoms.h"
#include "ForcesEnergiesCutoff.h"
#include "LatticeGen.h"
#include "berendsen.h"
#include "domain.h"
#include "gupta.h"
#include "metrics.h"
#include "mpi.h"
#include "neighbors.h"
#include "verlet.h"
#include "xyz.h"
#include <Eigen/Dense>
#include <iostream>

void write_values(std::ofstream &file, double kinetic_energy,
                  double potential_energy, double temperature, double strain,
                  double stress) {

    double energy = kinetic_energy + potential_energy;

    file << kinetic_energy << std::setw(15) << potential_energy << std::setw(15)
         << energy << std::setw(15) << temperature << std::setw(15) << strain
         << std::setw(15) << stress << std::endl;
}

double kinetic_energy_parallel(Atoms &atoms, Domain &domain) {
    double kinetic_energy = 0.0;
    for (int i = 1; i < domain.nb_local(); i++) {
        kinetic_energy += 0.5 *
                          pow(atoms.velocities.col(i).matrix().norm(), 2) *
                          atoms.masses(i);
    }
    return kinetic_energy;
}

double temperature_parallel(Atoms &atoms, Domain &domain) {
    return kinetic_energy_parallel(atoms, domain) /
           ((3.0 / 2.0) * BOLTZMANN * atoms.nb_atoms());
}

double avg_temp(double kin, int nb_atoms) {
    return kin / ((3.0 / 2.0) * BOLTZMANN * nb_atoms);
}

void berendsen_thermostat_parallel(Atoms &atoms, Domain &domain,
                                   double goal_temp, double timestep,
                                   double tau) {
    double lambda = std::sqrt(
        1.0 + (goal_temp / temperature_parallel(atoms, domain) - 1.0) *
                  timestep / tau);
    atoms.velocities *= lambda;
}

double compute_stress(Atoms &atoms, Domain domain, vec pots,
                      double border_width) {

    // borders of current sub domain
    auto border_left =
        domain.domain_length()[2] / domain.size() * domain.coordinate()[2];
    auto border_right =
        domain.domain_length()[2] / domain.size() * domain.coordinate()[2] +
        (domain.domain_length()[2] / domain.size());

    // borders plus minus the border_width
    auto border_left_left = border_left - border_width;
    auto border_right_right = border_right + border_width;

    double left_energies = 0.0;
    for (int i = 0; i < pots.rows(); i++) {
        if (atoms.positions(2, i) <= border_left &&
            atoms.positions(2, i) >= border_left_left) {
            left_energies += pots[i];
        }
    }
    double right_energies = 0.0;
    for (int i = 0; i < pots.rows(); i++) {
        if (atoms.positions(2, i) >= border_right &&
            atoms.positions(2, i) <= border_right_right) {
            right_energies += pots[i];
        }
    }
    std::cout << domain.rank() << " " << left_energies << " " << right_energies
              << std::endl;
    return std::min(left_energies, right_energies);
}

int main(int argc, char *argv[]) {

    std::cout << "Hello milestone 09!" << std::endl;
    MPI_Init(&argc, &argv);

    std::string path = "/home/sauerbac/MolecularDynamics/milestones/09/static/"
                       "whisker_small.xyz";
    std::string xyz_out =
        "/home/sauerbac/MolecularDynamics/milestones/09/static/xyz_out.xyz";
    std::ofstream file(xyz_out);
    std::string metrics = "/home/sauerbac/MolecularDynamics/milestones/09/"
                          "static/metrics.csv";
    std::ofstream metrics_file(metrics);

    Atoms system = read_atoms_no_velocities(path);
    Domain domain(MPI_COMM_WORLD, {50.0, 50.0, 144.1}, {1, 1, 6}, {0, 0, 1});
    double cutoff_distance = 4.0;
    NeighborList neighbor_list(cutoff_distance);

    double timestep = 1;
    double goal_temp = 300;
    int relax = 100;
    int period_steps = relax / timestep;

    double strain_zero = domain.domain_length()[2];
    double strain_add = 1.01;

    domain.enable(system);

    for (int i = 0; i < period_steps; i++) {
        verlet1(system, timestep);
        domain.exchange_atoms(system);
        domain.update_ghosts(system, 2 * cutoff_distance);
        neighbor_list.update(system);
        gupta_parallel(system, neighbor_list);
        verlet2(system, timestep);
        berendsen_thermostat_parallel(system, domain, goal_temp, timestep,
                                      relax / 100);
    }
    // start simu loop
    for (int j = 0; j < 100; j++) {

        // input energy
        domain.scale(
            system,
            Eigen::Array3d(50.0, 50.0, domain.domain_length()[2] * strain_add));

        // equilibrate
        for (int i = 0; i < period_steps; i++) {
            verlet1(system, timestep);
            domain.exchange_atoms(system);
            domain.update_ghosts(system, 2 * cutoff_distance);
            neighbor_list.update(system);
            gupta_parallel(system, neighbor_list);
            verlet2(system, timestep);
            berendsen_thermostat_parallel(system, domain, goal_temp, timestep,
                                          relax / 100);
        }

        // average
        vec potential_energies(period_steps);
        vec kinetic_energies(period_steps);
        vec stresses(period_steps);
        for (int i = 0; i < period_steps; i++) {
            verlet1(system, timestep);
            domain.exchange_atoms(system);
            domain.update_ghosts(system, 2 * cutoff_distance);

            neighbor_list.update(system);

            vec pots = gupta_parallel(system, neighbor_list);
            verlet2(system, timestep);

            // global potential energy
            double pot_global_temp;
            double pot = 0.0;
            for (int p = 0; p < domain.nb_local(); p++) {
                pot += pots.coeff(p);
            }
            MPI_Reduce(&pot, &pot_global_temp, 1, MPI_DOUBLE, MPI_SUM, 0,
                       MPI_COMM_WORLD);
            potential_energies.row(i) = pot_global_temp;

            // global kinteic energy
            double kin_global_temp;
            double kin = kinetic_energy_parallel(system, domain);
            MPI_Reduce(&kin, &kin_global_temp, 1, MPI_DOUBLE, MPI_SUM, 0,
                       MPI_COMM_WORLD);
            kinetic_energies.row(i) = kin_global_temp;

            // compute the local and gloabl stress
            double stress_global_temp;
            double stress = compute_stress(system, domain, potential_energies,
                                           2 * cutoff_distance);
            MPI_Reduce(&stress, &stress_global_temp, 1, MPI_DOUBLE, MPI_MIN, 0,
                       MPI_COMM_WORLD);
            stresses.row(i) = stress_global_temp;
        }

        // disable system to write to file
        domain.disable(system);
        if (domain.rank() == 0) {
            std::cout << j << std::endl;
            write_xyz(file, system);
            write_values(metrics_file, kinetic_energies.mean(),
                         potential_energies.mean(),
                         avg_temp(kinetic_energies.mean(), system.nb_atoms()),
                         (domain.domain_length()[2] - strain_zero) /
                             strain_zero,
                         stresses.mean());
        }
        domain.enable(system);
        domain.exchange_atoms(system);
        domain.update_ghosts(system, 2 * cutoff_distance);
    }

    MPI_Finalize();
    return 0;
}
