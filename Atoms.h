#ifndef __ATOMS_H
#define __ATOMS_H

#include "Types.h"
#include <Eigen/Dense>

class Atoms {
  public:
    mat positions;
    mat velocities;
    mat forces;
    vec masses;

    Atoms(const mat &p, const vec &m);
    Atoms(const mat &p, const mat &v, const vec &m);

    size_t nb_atoms() const;

    void verlet1(Atoms &atoms, double timestep);
    void verlet2(Atoms &atoms, double timestep);
};

#endif // __ATOMS_H