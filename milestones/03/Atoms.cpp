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
    assert_equal_length();
    return positions.cols();
}

void Atoms::resize(int len) {
    assert_equal_length();
    positions.conservativeResize(Eigen::NoChange_t::NoChange, len);
    velocities.conservativeResize(Eigen::NoChange_t::NoChange, len);
    forces.conservativeResize(Eigen::NoChange_t::NoChange, len);
    masses.conservativeResize(len, Eigen::NoChange_t::NoChange);
}

void Atoms::assert_equal_length() const {
    bool same_size = positions.cols() == velocities.cols() &&
                     velocities.cols() == forces.cols() &&
                     forces.cols() == masses.rows();
    if (!same_size) {
        std::cout << "positions:" << positions.cols() << std::endl;
        std::cout << "velocities:" << velocities.cols() << std::endl;
        std::cout << "forces:" << forces.cols() << std::endl;
        std::cout << "masses:" << masses.rows() << std::endl;
        std::cout << "1 is true, 0 is false:" << same_size << std::endl;
        throw std::runtime_error("Not same length in all arrays!");
    }
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
