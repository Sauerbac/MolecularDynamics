#include "verlet.h"
#include "Atoms.h"

void verlet1(Atoms &atoms, double timestep) {
    // TODO: check if multiplication is correct
    // atoms.velocities += 0.5 * atoms.forces * timestep / atoms.masses;
    atoms.velocities +=
        atoms.forces.rowwise() / atoms.masses.transpose() * 0.5 * timestep;
    atoms.positions += atoms.velocities * timestep;
}

void verlet2(Atoms &atoms, double timestep) {
    // TODO: check if multiplication is correct
    atoms.velocities +=
        atoms.forces.rowwise() / atoms.masses.transpose() * 0.5 * timestep;
}