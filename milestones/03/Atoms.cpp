#include "Atoms.h"

Atoms::Atoms(const mat &p, const vec &m)
    : positions{p},
      velocities{3, p.cols()},
      forces{3, p.cols()},
      masses{m} {
    velocities.setZero();
    forces.setZero();
}

Atoms::Atoms(const mat &p, const mat &v, const vec &m)
    : positions{p},
      velocities{v},
      forces{3, p.cols()},
      masses{m} {
    assert(p.cols() == v.cols());
    forces.setZero();
}

size_t Atoms::nb_atoms() const {
    return positions.cols();
}

void Atoms::verlet1(Atoms &atoms, double timestep) {
    // TODO: check if multiplication is correct
    atoms.velocities += 0.5 * atoms.forces * timestep / atoms.masses;
    atoms.positions += atoms.velocities * timestep;
}

void Atoms::verlet2(Atoms &atoms, double timestep) {
    // TODO: check if multiplication is correct
    atoms.velocities += 0.5 * atoms.velocities * timestep / atoms.masses;
}
