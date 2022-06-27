#include "Atoms.h"
#include "ForcesEnergies.h"
#include "neighbors.h"
#include <iostream>

int main() {
    std::cout << "Hello milestone 06!" << std::endl;

    str_vec names{{"H", "H", "H", "H"}};
    mat positions(3, 4);
    positions << 0, 7, 0, 0, 0, 0, 7, 6, 0, 0, 0, 0;

    Atoms atoms(positions, names);
    NeighborList NeighborList(0.5);

    return 0;
}
