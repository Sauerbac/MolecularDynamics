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

int main() {
    std::cout << "Hello milestone 07!" << std::endl;

    std::string path =
        "/home/sauerbach/HPC/MolecularDynamics/static/cluster_309.xyz";
    Atoms system = read_atoms_no_velocities(path);
    // Atoms system = cubic(300, 3.0, 196.966569, 0.01);
    std::cout << system.nb_atoms() << std::endl;

    // std::string xyz_out =
    //     "/home/sauerbach/HPC/MolecularDynamics/static/out_no_domain.xyz";
    // std::ofstream file(xyz_out);

    // std::string metrics =
    //     "/home/sauerbach/HPC/MolecularDynamics/static/metrics_309.csv";
    // std::ofstream metrics_file(metrics);

    std::string avg_temp_path =
        "/home/sauerbach/HPC/MolecularDynamics/static/temps_309.csv";
    std::ofstream temp_file(avg_temp_path);

    NeighborList neighbor_list(5.0);

    double timestep = 1;
    double period_time = 10000;
    int period_steps = period_time / timestep;
    int tau = 100;

    double temp = 460;

    for (int i = 0; i < period_steps; i++) {

        neighbor_list.update(system);
        verlet1(system, timestep);
        double pot = gupta(system, neighbor_list);
        verlet2(system, timestep);
        berendsen_thermostat(system, temp, timestep, tau);
    }

    for (int j = 0; j < 100; j++) {

        vec kinetic_energys(period_steps);
        double pot = 0;

        for (int i = 0; i < period_steps; i++) {

            neighbor_list.update(system);
            verlet1(system, timestep);
            pot = gupta(system, neighbor_list);
            verlet2(system, timestep);

            kinetic_energys.row(i) = kinetic_energy(system);
        }

        double a_temp = avg_temp(kinetic_energys.mean(), system.nb_atoms());
        write_avgs(temp_file, a_temp, pot);
        std::cout << "Sclaing" << std::endl;
        double factor = sqrt(1 / kinetic_energys.mean() + 1);
        system.velocities *= factor;

        for (int i = 0; i < period_steps; i++) {

            neighbor_list.update(system);
            verlet1(system, timestep);
            pot = gupta(system, neighbor_list);
            verlet2(system, timestep);

            // write_metrics(metrics_file, system, pot);
            // if (i % 100 == 0) {
            //     write_xyz(file, system);
            // }
        }
    }
    std::cout << period_steps << std::endl;

    // file.close();
}
