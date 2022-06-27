#include "ForcesEnergiesCutoff.h"

void lj_neighbors(Atoms &atoms, NeighborList &neighbor_list, double epsilon,
                  double sigma) {
    atoms.forces.setZero();
    for (auto [i, j] : neighbor_list) {
        if (i != j) {

            double distance = (atoms.positions.col(i) - atoms.positions.col(j))
                                  .matrix()
                                  .norm();

            // Compute the new forces
            double pauli_force = 12.0 * pow(sigma, 12.0) / pow(distance, 13.0);
            double london_force = 6.0 * pow(sigma, 6.0) / pow(distance, 7.0);
            Eigen::Vector3d single_force =
                4.0 * epsilon * (pauli_force - london_force) *
                (atoms.positions.col(i) - atoms.positions.col(j))
                    .matrix()
                    .normalized();
            atoms.forces.col(i) += single_force.array();
        }
    }
}

double potential_energy_cutoff(Atoms &atoms, NeighborList &neighbor_list,
                               double epsilon, double sigma) {
    double cutoff = neighbor_list.get_interaction_range();
    double pauli_energy = pow(sigma / cutoff, 12.0);
    double london_energy = pow(sigma / cutoff, 6.0);
    double cutoff_potential = 4.0 * epsilon * (pauli_energy - london_energy);

    double pot_energy = 0;
    for (auto [i, j] : neighbor_list) {
        if (i < j) {
            double distance = (atoms.positions.col(i) - atoms.positions.col(j))
                                  .matrix()
                                  .norm();
            // Compute potential Energy
            double pauli_energy = pow(sigma / distance, 12.0);
            double london_energy = pow(sigma / distance, 6.0);
            pot_energy += 4.0 * epsilon * (pauli_energy - london_energy) -
                          cutoff_potential;
        }
    }
    return 0.5 * pot_energy;
}
