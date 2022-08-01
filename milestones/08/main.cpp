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

    std::string path =
        "/home/sauerbach/HPC/MolecularDynamics/static/cluster_923.xyz";
    std::string xyz_out =
        "/home/sauerbach/HPC/MolecularDynamics/static/out2.xyz";
    std::ofstream file(xyz_out);
    std::string metrics =
        "/home/sauerbach/HPC/MolecularDynamics/static/metrics.csv";
    std::ofstream metrics_file(metrics);

    Atoms system = read_atoms_no_velocities(path);
    Domain domain(MPI_COMM_WORLD, {30.0, 30.0, 30.0}, {1, 1, 4}, {0, 0, 0});
    double cutoff_distance = 10.0;
    NeighborList neighbor_list(cutoff_distance);

    double timestep = 10;
    double goal_temp = 500;
    double relaxation_time = 10000;

    int num_steps = relaxation_time / timestep;
    domain.enable(system);
    domain.exchange_atoms(system);
    domain.update_ghosts(system, 2 * cutoff_distance);
    for (int i = 0; i < num_steps; i++) {
        // std::cout << domain.rank() << "  " << system.positions.cols() << "  "
        //           << system.nb_atoms() << std::endl;
        if (system.nb_atoms() > 0) {
            neighbor_list.update(system);
        }
        verlet1(system, timestep);
        domain.exchange_atoms(system);
        domain.update_ghosts(system, 2 * cutoff_distance);
        std::cout << i << "   " << domain.rank() << std::endl;

        double pot = 0.0;
        if (system.nb_atoms() > 0) {
            pot = gupta(system, neighbor_list);
        }
        verlet2(system, timestep);
        // berendsen_thermostat(system, goal_temp, timestep,
        // relaxation_time); std::cout << i << std::endl;
        if (i % num_steps / 100 == 0) {
            domain.disable(system);
            write_xyz(file, system);
            write_metrics(metrics_file, system, pot);
            domain.enable(system);
            domain.exchange_atoms(system);
            domain.update_ghosts(system, 2 * cutoff_distance);
        }
    }

    MPI_Finalize();
    return 0;
}
