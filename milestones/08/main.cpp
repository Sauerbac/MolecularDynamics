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

int main(int argc, char *argv[]) {

    std::cout << "Hello milestone 08!" << std::endl;
    MPI_Init(&argc, &argv);

    std::string path = "/home/sauerbac/MolecularDynamics/milestones/08/static/"
                       "cluster_3871.xyz";
    std::string xyz_out =
        "/home/sauerbac/MolecularDynamics/milestones/08/static/out_8.xyz";
    std::ofstream file(xyz_out);
    std::string metrics = "/home/sauerbac/MolecularDynamics/milestones/08/"
                          "static/metrics_122.csv";
    std::ofstream metrics_file(metrics);

    Atoms system = read_atoms_no_velocities(path);
    Domain domain(MPI_COMM_WORLD, {50.0, 50.0, 50.0}, {1, 2, 2}, {0, 0, 0});
    double cutoff_distance = 4.0;
    NeighborList neighbor_list(cutoff_distance);

    double timestep = 1;
    double goal_temp = 1000;
    double relaxation_time = 50000;

    int num_steps = relaxation_time / timestep;

    domain.enable(system);
    domain.exchange_atoms(system);
    domain.update_ghosts(system, 2 * cutoff_distance);
    for (int i = 0; i < num_steps; i++) {
        verlet1(system, timestep);
        domain.exchange_atoms(system);
        domain.update_ghosts(system, 2 * cutoff_distance);

        // std::cout << "Domain " << domain.rank()
        //           << ", number atoms: " << system.nb_atoms() << std::endl;

        double pot = 0.0;
        if (system.nb_atoms() > 0) {
            neighbor_list.update(system);

            vec potential_energies = gupta_parallel(system, neighbor_list);
            for (int p = 0; p < domain.nb_local(); p++) {
                pot += potential_energies.coeff(p);
            }
            verlet2(system, timestep);
            // berendsen_thermostat(system, goal_temp, timestep,
            //                      relaxation_time / 100);
        }
        double pot_global;
        MPI_Reduce(&pot, &pot_global, 1, MPI_DOUBLE, MPI_SUM, 0,
                   MPI_COMM_WORLD);
        if (i % (num_steps / 1000) == 0) {

            domain.disable(system);
            if (domain.rank() == 0) {
                std::cout << system.nb_atoms() << std::endl;
                write_xyz(file, system);
                write_metrics(metrics_file, system, pot_global);
            }
            domain.enable(system);
            domain.exchange_atoms(system);
            domain.update_ghosts(system, 2 * cutoff_distance);
        }
    }

    MPI_Finalize();
    return 0;
}
