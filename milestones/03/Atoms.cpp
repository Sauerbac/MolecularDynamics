#include "Atoms.h"
#include "constants.h"
#include <iostream>

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

Atoms::Atoms(int nb_atoms)
    : positions{3, nb_atoms},
      velocities{3, nb_atoms},
      forces{3, nb_atoms},
      masses{nb_atoms} {
    positions.setZero();
    velocities.setZero();
    forces.setZero();
}

Atoms::Atoms(mat &p, str_vec &n) {
    positions = p;
    velocities = mat::Zero(3, p.cols());
    forces = mat::Zero(3, p.cols());
    masses = name2masses(n);
}

Atoms::Atoms(mat &p, const mat &v, str_vec &n) {
    positions = p;
    velocities = v;
    forces = mat::Zero(3, p.cols());
    masses = name2masses(n);
}

size_t Atoms::nb_atoms() const {
    return positions.cols();
}

vec name2masses(str_vec &names) {
    int len = names.rows();

    vec masses{len};
    for (int i = 0; i < len; i++) {
        std::string atom_name = names[i];
        double mass = ATOMIC_MASSES.at(atom_name);

        masses.row(i) = mass * MASSFACTOR;
    }
    return masses;
}
