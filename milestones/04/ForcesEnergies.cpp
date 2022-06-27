#include "ForcesEnergies.h"

void lj_direct_summation(Atoms &atoms, double epsilon, double sigma) {
    for (int i = 0; i < atoms.nb_atoms(); i++) {
        Eigen::Vector3d force_sum = Eigen::Vector3d::Zero();
        for (int j = 0; j < atoms.nb_atoms(); j++) {
            if (i != j) {
                double distance =
                    (atoms.positions.col(i) - atoms.positions.col(j))
                        .matrix()
                        .norm();

                // Compute the new forces
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
            }
        }
        atoms.forces.col(i) = force_sum;
    }
}

double kinetic_energy(Atoms &atoms) {
    double kinetic_energy = 0.0;
    for (int i = 1; i < atoms.nb_atoms(); i++) {
        kinetic_energy += 0.5 *
                          pow(atoms.velocities.col(i).matrix().norm(), 2) *
                          atoms.masses(i);
    }
    return kinetic_energy;
}

double potential_energy(Atoms &atoms, double epsilon, double sigma) {
    double pot_energy = 0;
    for (int i = 0; i < atoms.nb_atoms(); i++) {
        for (int j = 0; j < atoms.nb_atoms(); j++) {
            if (i != j) {
                double distance =
                    (atoms.positions.col(i) - atoms.positions.col(j))
                        .matrix()
                        .norm();
                // Compute potential Energy
                double pauli_energy = pow(sigma / distance, 12.0);
                double london_energy = pow(sigma / distance, 6.0);
                pot_energy += 4.0 * epsilon * (pauli_energy - london_energy);
            }
        }
    }
    return 0.5 * pot_energy;
}
double total_energy(Atoms &atoms, double sigma, double epsilon) {
    return potential_energy(atoms, sigma, epsilon) + kinetic_energy(atoms);
}
double temperature(Atoms &atoms) {
    return kinetic_energy(atoms) / ((3.0 / 2.0) * BOLTZMANN * atoms.nb_atoms());
}
