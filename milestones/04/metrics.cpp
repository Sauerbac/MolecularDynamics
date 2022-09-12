#include "metrics.h"

std::tuple<double, double, double, double>
calc_metrics(Atoms &atoms, double potential_energy) {
    double kinetic = kinetic_energy(atoms);
    double energy = kinetic + potential_energy;
    double temp = temperature(atoms);
    return std::tuple(potential_energy, kinetic, energy, temp);
}

void write_metrics(std::ofstream &file, Atoms &atoms, double potential_energy) {

    double kinetic = kinetic_energy(atoms);
    double energy = kinetic + potential_energy;
    double temp = temperature(atoms);

    file << kinetic << std::setw(15) << potential_energy << std::setw(15)
         << energy << std::setw(15) << temp << std::endl;
}
