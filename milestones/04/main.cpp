#include "Atoms.h"
#include "ForcesEnergies.h"
#include "metrics.h"
#include "verlet.h"
#include "xyz.h"
#include <iostream>

int main() {
    std::cout << "Hello milestone 05!" << std::endl;

    std::string path = "/home/sauerbach/HPC/MolecularDynamics/static/lj54.xyz";
    Atoms system = read_atoms_no_velocities(path);

    std::string output =
        "/home/sauerbach/HPC/MolecularDynamics/milestones/04/static/0.45.xyz";
    std::ofstream file(output);
    std::string metrics =
        "/home/sauerbach/HPC/MolecularDynamics/milestones/04/static/0.45.csv";
    std::ofstream metrics_file(metrics);

    double epsilon = 1.0;
    double sigma = 1.0;
    double mass = 1.0;
    int simu_length = 1000;
    double timestep = 0.45;
    int steps = simu_length / timestep;
    int num_snapshots = 1000;

    std::cout << simu_length << std::endl;
    std::cout << timestep << std::endl;
    std::cout << steps << std::endl;
    std::cout << num_snapshots << std::endl;

    for (int i = 0; i < steps; i++) {
        verlet1(system, timestep);
        lj_direct_summation(system, epsilon, sigma);
        verlet2(system, timestep);
        // std::cout << temperature(system) << std::endl;

        if (i % (steps / num_snapshots) == 0) {

            std::cout << i * timestep << std::endl;
            write_xyz(file, system);
            double pot = potential_energy(system, sigma, epsilon);
            write_metrics(metrics_file, system, pot);
            // std::cout << system.forces << std::endl;
        }
    }

    file.close();

    return 0;
}
