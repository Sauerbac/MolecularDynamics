#ifndef __ATOMS_H
#define __ATOMS_H

#include <Eigen/Dense>

class Atoms {
  public:
    Eigen::Array3Xd positions;
    Eigen::Array3Xd velocities;
    Eigen::Array3Xd forces;
    Eigen::ArrayXd masses;

    Atoms(const Eigen::Array3Xd &p, const Eigen::ArrayXd &m);
    Atoms(const Eigen::Array3Xd &p, const Eigen::Array3Xd &v,
          const Eigen::ArrayXd &m);

    size_t nb_atoms() const;

    void verlet1(Atoms &atoms, double timestep);
    void verlet2(Atoms &atoms, double timestep);
};

#endif // __ATOMS_H