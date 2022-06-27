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
    Atoms(int nb_atoms);
    Atoms(mat &p, str_vec &n);

    size_t nb_atoms() const;
};

vec name2masses(str_vec &names);

#endif // __ATOMS_H