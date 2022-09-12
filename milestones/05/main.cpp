// #include "Atoms.h"
// #include "ForcesEnergies.h"
// #include "LatticeGen.h"
// #include "berendsen.h"
// #include "metrics.h"
// #include "verlet.h"
// #include "xyz.h"
// #include <Eigen/Dense>
// #include <iostream>

// int main() {
//     std::cout << "Hello milestone 05!" << std::endl;

//     std::string path = "/home/sauerbach/HPC/MolecularDynamics/milestones/05/"
//                        "static/example.xyz";
//     // Atoms system = cubic(50, 1.0, 1.008, 0.1);
//     Atoms system = read_atoms(path);

//     std::string output =
//         "/home/sauerbach/HPC/MolecularDynamics/milestones/05/static/test.xyz";
//     std::ofstream file(output);
//     std::string metrics =
//         "/home/sauerbach/HPC/MolecularDynamics/milestones/05/static/test.csv";
//     std::ofstream metrics_file(metrics);

//     double epsilon = 1.0;
//     double sigma = 1.0;
//     double mass = 1.0;
//     int simu_length = 10;
//     double timestep = 0.01;
//     int steps = simu_length / timestep;
//     int num_snapshots = 100;

//     double goal_temp = 0.2;
//     double tau = 1;

//     system.masses = 1.008;

//     std::cout << temperature(system) << std::endl;
//     std::cout << system.velocities << std::endl;

//     for (int i = 0; i < steps; i++) {
//         // if (i < steps * 0.1) {

//         //     verlet1(system, timestep);
//         //     lj_direct_summation(system, epsilon, sigma);
//         //     verlet2(system, timestep);
//         // }

//         std::cout << temperature(system) << std::endl;
//         berendsen_thermostat(system, goal_temp, timestep, tau);
//         std::cout << temperature(system) << std::endl;
//         if (i % (steps / num_snapshots) == 0) {

//             std::cout << i * timestep << std::endl;
//             std::cout << temperature(system) << std::endl;
//             write_xyz(file, system);
//             double pot = potential_energy(system, sigma, epsilon);
//             write_metrics(metrics_file, system, pot);
//             // std::cout << system.forces << std::endl;
//         }
//     }

//     file.close();

//     return 0;
// }

// used to measure time complexity
#include "Atoms.h"
#include "ForcesEnergies.h"
#include "LatticeGen.h"
#include "berendsen.h"
#include "metrics.h"
#include "verlet.h"
#include "xyz.h"
#include <Eigen/Dense>
#include <chrono>
#include <iostream>
using namespace std::chrono;

int main(int argc, char **argv) {
    // std::cout << "Hello milestone 05!" << std::endl;
    // std::cout << "num atoms:" << std::endl;
    int num_atoms = std::stoi(argv[1]);
    // std::cout << num_atoms << std::endl;

    Atoms system = cubic(num_atoms, 1.1, 1.008, 0.1);
    auto start = high_resolution_clock::now();

    double timestep = 0.01;
    double epsilon = 1.0;
    double sigma = 1.0;

    for (int i = 0; i < 10000; i++) {
        verlet1(system, timestep);
        lj_direct_summation(system, epsilon, sigma);
        verlet2(system, timestep);
    }
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    std::cout << duration.count() << std::endl;

    return 0;
}