// #include "Atoms.h"
// #include "ForcesEnergiesCutoff.h"
// #include "LatticeGen.h"
// #include "berendsen.h"
// #include "neighbors.h"
// #include "verlet.h"
// #include "xyz.h"
// #include <Eigen/Dense>
// #include <iostream>

// int main() {
//     std::cout << "Hello milestone 05!" << std::endl;

//     std::string path =
//     "/home/sauerbach/HPC/MolecularDynamics/static/lj54.xyz";
//     // Atoms system = read_atoms(path);
//     Atoms system = cubic(1000, 1.0, 1.008, 0.01);
//     std::cout << system.nb_atoms() << std::endl;

//     std::string output =
//     "/home/sauerbach/HPC/MolecularDynamics/static/out.xyz"; std::ofstream
//     file(output);

//     double timestep = 0.01;
//     double epsilon = 1.0;
//     double sigma = 1.0;
//     double goal_temp = 1000.0;
//     double relaxation_time = 1.0;

//     NeighborList neighbor_list(2.0);

//     for (int i = 0; i < 100; i++) {
//         write_xyz(file, system);

//         neighbor_list.update(system);
//         verlet1(system, timestep);
//         lj_neighbors(system, neighbor_list, 1.0, 1.0);
//         verlet2(system, timestep);
//         berendsen_thermostat(system, goal_temp, timestep, relaxation_time);
//         std::cout << temperature(system) << std::endl;
//     }

//     file.close();

//     return 0;
// }

#include "Atoms.h"
#include "ForcesEnergies.h"
#include "ForcesEnergiesCutoff.h"
#include "LatticeGen.h"
#include "berendsen.h"
#include "metrics.h"
#include "neighbors.h"
#include "verlet.h"
#include "xyz.h"
#include <Eigen/Dense>
#include <chrono>
#include <iostream>
using namespace std::chrono;

int main(int argc, char **argv) {
    // std::cout << "Hello milestone 05!" << std::endl;
    // std::cout << "num atoms:" << std::endl;
    // int num_atoms = 500;
    int num_atoms = std::stoi(argv[1]);
    // std::cout << num_atoms << std::endl;

    // std::string output =
    // "/home/sauerbach/HPC/MolecularDynamics/static/out.xyz"; std::ofstream
    // file(output); std::string metrics =
    //     "/home/sauerbach/HPC/MolecularDynamics/static/metrics.csv";
    // std::ofstream metrics_file(metrics);

    Atoms system = cubic(num_atoms, 1.1, 1.008, 0.1);
    auto start = high_resolution_clock::now();

    double timestep = 0.01;
    double epsilon = 1.0;
    double sigma = 1.0;
    int steps = 10000;

    NeighborList neighbor_list(3.0);

    for (int i = 0; i < steps; i++) {
        neighbor_list.update(system);
        verlet1(system, timestep);
        lj_neighbors(system, neighbor_list, epsilon, sigma);
        verlet2(system, timestep);
        // if (i % (steps / 100) == 0) {

        //     std::cout << i * timestep << std::endl;
        //     write_xyz(file, system);
        //     double pot = potential_energy(system, sigma, epsilon);
        //     write_metrics(metrics_file, system, pot);
        //     // std::cout << system.forces << std::endl;
        // }
    }
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    std::cout << duration.count() << std::endl;

    return 0;
}
