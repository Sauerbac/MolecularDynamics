#include "LatticeGen.h"
#include <iostream>

Atoms cubic(int num_atoms, double side_length, double mass_all_atoms,
            double max_offset) {
    int num_atoms_side = std::ceil(std::cbrt(num_atoms)) + 2;
    int max_atoms = std::pow(num_atoms_side, 3);

    mat positions{3, max_atoms};
    for (int i = 0; i < num_atoms_side; i++) {
        for (int j = 0; j < num_atoms_side; j++) {
            for (int k = 0; k < num_atoms_side; k++) {
                int index = i * num_atoms_side * num_atoms_side +
                            j * num_atoms_side + k;

                positions.col(index)
                    << i * side_length + random_between_0_and(max_offset),
                    j * side_length + random_between_0_and(max_offset),
                    k * side_length + random_between_0_and(max_offset);
            }
        }
    }

    vec center_of_mass = positions.rowwise().mean();
    for (int i = 0; i < max_atoms - num_atoms; i++) {
        int furthest_index = 0;
        double furthest_distance = 0.0;
        for (int j = 0; j < positions.cols(); j++) {
            double distance =
                (positions.col(j) - center_of_mass).matrix().norm();
            if (distance > furthest_distance) {
                furthest_distance = distance;
                furthest_index = j;
            }
        }
        removeColumn(positions, furthest_index);
    }

    vec masses = vec::Constant(num_atoms, mass_all_atoms);
    return Atoms(positions, masses);
}

// https://stackoverflow.com/questions/686353/random-float-number-generation
double random_between_0_and(double high) {
    return static_cast<double>(rand()) / (static_cast<float>(RAND_MAX / high));
}

// https://stackoverflow.com/questions/13290395/how-to-remove-a-certain-row-or-column-while-using-eigen-library-c
void removeColumn(mat &matrix, unsigned int colToRemove) {
    unsigned int numRows = matrix.rows();
    unsigned int numCols = matrix.cols() - 1;

    if (colToRemove < numCols)
        matrix.block(0, colToRemove, numRows, numCols - colToRemove) =
            matrix.block(0, colToRemove + 1, numRows, numCols - colToRemove);

    matrix.conservativeResize(numRows, numCols);
}