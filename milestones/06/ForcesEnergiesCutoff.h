#ifndef __FORCESENERGIES_NEIGHBORS_H
#define __FORCESENERGIES_NEIGHBORS_H

#include "Atoms.h"
#include "constants.h"
#include "neighbors.h"
#include <Eigen/Dense>

void lj_neighbors(Atoms &atoms, NeighborList &neighbor_list,
                  double epsilon = 1.0, double sigma = 1.0);

double potential_energy_cutoff(Atoms &atoms, NeighborList &neighbor_list,
                               double epsilon, double sigma);

#endif // __FORCESENERGIES_NEIGHBORS_H
