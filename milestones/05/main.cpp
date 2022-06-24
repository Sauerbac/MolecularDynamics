#include "Atoms.h"
#include "ForcesEnergies.h"
#include "LatticeGen.h"
#include "berendsen.h"
#include "verlet.h"
#include "xyz.h"
#include <Eigen/Dense>
#include <iostream>

int main() {
    std::cout << "Hello milestone 05!" << std::endl;

    std::string path = "/home/sauerbach/HPC/MolecularDynamics/static/lj54.xyz";
    // Atoms system = read_atoms(path);
    Atoms system = cubic(300, 1.0, 1.008, 0.01);
    std::cout << system.nb_atoms() << std::endl;

    std::string output = "/home/sauerbach/HPC/MolecularDynamics/static/out.xyz";
    std::ofstream file(output);

    double timestep = 0.01;
    double epsilon = 1.0;
    double sigma = 1.0;
    double goal_temp = 1000.0;
    double relaxation_time = 1.0;

    for (int i = 0; i < 1000; i++) {
        write_xyz(file, system);
        verlet1(system, timestep);
        lj_direct_summation(system, 1.0, 1.0);
        verlet2(system, timestep);
        berendsen_thermostat(system, goal_temp, timestep, relaxation_time);
        std::cout << temperature(system) << std::endl;
    }

    file.close();

    return 0;
}
