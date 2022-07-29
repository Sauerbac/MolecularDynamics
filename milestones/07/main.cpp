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

int main() {
    std::cout << "Hello milestone 07!" << std::endl;

    std::string path =
        "/home/sauerbach/HPC/MolecularDynamics/static/cluster_923.xyz";
    Atoms system = read_atoms_no_velocities(path);
    // Atoms system = cubic(300, 3.0, 196.966569, 0.01);
    std::cout << system.nb_atoms() << std::endl;

    std::string xyz_out =
        "/home/sauerbach/HPC/MolecularDynamics/static/out.xyz";
    std::ofstream file(xyz_out);

    std::string metrics =
        "/home/sauerbach/HPC/MolecularDynamics/static/metrics.csv";
    std::ofstream metrics_file(metrics);

    NeighborList neighbor_list(7.0);

    double timestep = 10;
    double goal_temp = 500;
    double relaxation_time = 10000;

    int num_steps = relaxation_time / timestep;

    for (int j = 0; j < 50; j++) {
        for (int i = 0; i < num_steps; i++) {
            if (i % 10 == 0) {
                write_xyz(file, system);
            }

            neighbor_list.update(system);
            verlet1(system, timestep);
            double pot = gupta(system, neighbor_list);
            verlet2(system, timestep);
            berendsen_thermostat(system, goal_temp, timestep, relaxation_time);
            // std::cout << temperature(system) << std::endl;
            write_metrics(metrics_file, system, pot);
        }
        goal_temp = goal_temp + 10;
        for (int i = 0; i < num_steps * 2; i++) {
            if (i % 10 == 0) {
                write_xyz(file, system);
            }

            neighbor_list.update(system);
            verlet1(system, timestep);
            double pot = gupta(system, neighbor_list);
            verlet2(system, timestep);
            // std::cout << temperature(system) << std::endl;
            write_metrics(metrics_file, system, pot);
        }
        std::cout << j << std::endl;
        std::cout << goal_temp << ' ' << temperature(system) << std::endl;
    }

    file.close();

    return 0;
}
