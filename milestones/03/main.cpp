#include "Atoms.h"
#include <iostream>

int main() {
    std::cout << "Hello milestone 03!" << std::endl;

    Eigen::ArrayXd vec = Eigen::ArrayXd((1));
    Eigen::Array3Xd mat(3, 1);

    mat << 0, 0, 0;
    Atoms atom = Atoms(mat, vec);

    std::cout << atom.nb_atoms() << std::endl;

    return 0;
}
