#include "Atoms.h"
#include "ForcesEnergies.h"
#include <fstream>
#include <iomanip>
#include <iostream>

std::tuple<double, double, double, double>
calc_metrics(Atoms &atoms, double potential_energy);

void write_metrics(std::ofstream &file, Atoms &atoms, double potential_energy);
