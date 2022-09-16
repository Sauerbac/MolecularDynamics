#ifndef __FORCESENERGIES_H
#define __FORCESENERGIES_H

#include "Atoms.h"
#include "constants.h"
#include <Eigen/Dense>

void lj_direct_summation(Atoms &atoms, double epsilon = 1.0,
                         double sigma = 1.0);

double lj_direct_summation_test(Atoms &atoms, double epsilon = 1.0,
                                double sigma = 1.0);

double kinetic_energy(Atoms &atoms);
double potential_energy(Atoms &atoms, double sigma, double epsilon);
double total_energy(Atoms &atoms, double sigma, double epsilon);
double temperature_lj(Atoms &atoms);
double temperature(Atoms &atoms);

#endif // __FORCESENERGIES_H
