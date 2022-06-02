#ifndef __FORCESENERGIES_H
#define __FORCESENERGIES_H

#include "Atoms.h"

double lj_direct_summation(Atoms &atoms, double epsilon = 1.0,
                           double sigma = 1.0);

#endif // __FORCESENERGIES_H
