#include "berendsen.h"

void berendsen_thermostat(Atoms &atoms, double goal_temp, double timestep,
                          double tau) {
    double lambda = std::sqrt(1.0 + (goal_temp / temperature(atoms) - 1.0) *
                                        timestep / tau);
    atoms.velocities *= lambda;
}