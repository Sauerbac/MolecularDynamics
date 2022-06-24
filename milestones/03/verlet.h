#ifndef __VERLET_H
#define __VERLET_H

#include "Atoms.h"

void verlet1(Atoms &atoms, double timestep);
void verlet2(Atoms &atoms, double timestep);

#endif // __VERLET_H