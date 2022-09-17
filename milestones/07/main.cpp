#include "Atoms.h"
#include "ForcesEnergiesCutoff.h"
#include "LatticeGen.h"
#include "berendsen.h"
#include "gupta.h"
#include "metrics.h"
#include "neighbors.h"
#include "verlet.h"
#include "xyz.h"
#include <Eigen/Dense>
#include <iostream>

double avg_temp(double avg_kin, int nb_atoms) {
    return avg_kin / ((3.0 / 2.0) * BOLTZMANN * nb_atoms);
}

void write_avgs(std::ofstream &file, double avg_temp, double potential_energy) {
    file << avg_temp << std::setw(15) << potential_energy << std::endl;
}

double std_dev(vec vec) {
    return std::sqrt((vec - vec.mean()).square().sum() / (vec.size() - 1));
}

Atoms load_system() {
    std::string path = "/home/sauerbac/MolecularDynamics/milestones/07/static/"
                       "cluster_3871.xyz";
    Atoms system = read_atoms_no_velocities(path);
    // Atoms system = cubic(300, 3.0, 196.966569, 0.01);
    std::cout << system.nb_atoms() << std::endl;
    return system;
}

int main() {
    auto system = load_system();
    std::string output =
        "/home/sauerbac/MolecularDynamics/milestones/07/static/out.xyz";
    std::ofstream file(output);

    std::string metrics =
        "/home/sauerbac/MolecularDynamics/milestones/07/static/metrics.csv";
    std::ofstream metrics_file(metrics);

    std::string avg_temp_path = "/home/sauerbac/MolecularDynamics/milestones/"
                                "07/static/temps_3871.csv";
    std::ofstream temp_file(avg_temp_path);

    NeighborList neighbor_list(8.0);

    double timestep = 2;
    double period_time = 2000;
    int period_steps = period_time / timestep;

    double temp = 500;
    int tau = 100;

    // let system equibrilate
    for (int i = 0; i < 2 * period_steps; i++) {
        neighbor_list.update(system);
        verlet1(system, timestep);
        double pot = gupta(system, neighbor_list);
        verlet2(system, timestep);
        berendsen_thermostat(system, temp, timestep, tau);
    }
    // increment x times
    for (int j = 0; j < 200; j++) {
        vec kinetic_energys(period_steps);
        double pot = 0;

        // scale temp, add Delta Q
        double factor = sqrt((30 / kinetic_energy(system)) + 1);
        system.velocities *= factor;

        for (int i = 0; i < period_steps; i++) {
            neighbor_list.update(system);
            verlet1(system, timestep);
            pot = gupta(system, neighbor_list);
            verlet2(system, timestep);
            // metrics_file << temperature(system) << std::endl;
        }

        for (int i = 0; i < period_steps; i++) {
            neighbor_list.update(system);
            verlet1(system, timestep);
            pot = gupta(system, neighbor_list);
            verlet2(system, timestep);
            // metrics_file << temperature(system) << std::endl;

            kinetic_energys.row(i) = kinetic_energy(system);
        }

        double a_temp = avg_temp(kinetic_energys.mean(), system.nb_atoms());
        write_avgs(temp_file, a_temp, kinetic_energy(system) + pot);

        std::cout << "current kin: " << kinetic_energy(system)
                  << " avg_kin: " << kinetic_energys.mean()
                  << " std: " << std_dev(kinetic_energys) << std::endl;
        std::cout << "Run: " << j << " current temp: " << temperature(system)
                  << " avg_temp: " << a_temp << std::endl;
    }
    std::cout << period_steps << std::endl;

    metrics_file.close();
    temp_file.close();
    file.close();
}

// 923, 2, 2000, 10
// 561, 2, 2000, 6
// 309, 2, 2000, 3.3