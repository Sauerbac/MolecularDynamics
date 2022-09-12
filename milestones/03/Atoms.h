#ifndef __ATOMS_H
#define __ATOMS_H

#include "Types.h"
#include <Eigen/Dense>

struct Atoms {
    mat positions;
    mat velocities;
    mat forces;
    vec masses;

    Atoms(const mat &p, const vec &m);
    Atoms(const mat &p, const mat &v, const vec &m);
    Atoms(int nb_atoms);
    Atoms(mat &p, str_vec &n);
    Atoms(mat &p, const mat &v, str_vec &n);

    size_t nb_atoms() const;
    void resize(int len);
    void assert_equal_length() const;
};

vec name2masses(str_vec &names);
void removeRow(Eigen::MatrixXd &matrix, unsigned int rowToRemove);
void removeColumn(Eigen::MatrixXd &matrix, unsigned int colToRemove);

#endif // __ATOMS_H