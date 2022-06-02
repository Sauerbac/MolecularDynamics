#include "ForcesEnergies.h"
#include <Eigen/Dense>
#include <math.h>

double lj_direct_summation(Atoms &atoms, double epsilon, double sigma) {
    double pot_energy = 0;
    for (int i = 0; i < atoms.nb_atoms(); i++) {
        Eigen::Vector3d force_sum = Eigen::Vector3d::Zero();
        for (int j = 0; j < atoms.nb_atoms(); j++) {
            if (i != j) {
                double distance =
                    (atoms.positions.col(i) - atoms.positions.col(j))
                        .matrix()
                        .norm();

                // first, calculate the new forces
                double pauli_force =
                    12.0 * pow(sigma, 12.0) / pow(distance, 13.0);
                double london_force =
                    6.0 * pow(sigma, 6.0) / pow(distance, 7.0);
                Eigen::Vector3d single_force =
                    4.0 * epsilon * (pauli_force - london_force) *
                    (atoms.positions.col(i) - atoms.positions.col(j))
                        .matrix()
                        .normalized();
                force_sum += single_force;

                // second, compute the potential energy
                double pauli_energy = pow(sigma / distance, 12.0);
                double london_energy = pow(sigma / distance, 6.0);
                pot_energy += 4.0 * epsilon * (pauli_energy - london_energy);
            }
        }
        atoms.forces.col(i) = force_sum;
    }
    return 0.5 * pot_energy;
}
