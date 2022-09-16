#ifndef BERENDSEN_H
#define BERENDSEN_H

#include "Atoms.h"
#include "ForcesEnergies.h"

void berendsen_thermostat_lj(Atoms &atoms, double goal_temp, double timestep,
                             double tau);
void berendsen_thermostat(Atoms &atoms, double goal_temp, double timestep,
                          double tau);

#endif // BERENDSEN_H