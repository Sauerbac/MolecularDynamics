#include "Atoms.h"
#include <Eigen/Dense>

Atoms cubic(int num_atoms, double side_length, double mass_all_atoms,
            double offset_from_perfect_lattice);
double random_between_0_and(double high);
void removeColumn(mat &matrix, unsigned int colToRemove);