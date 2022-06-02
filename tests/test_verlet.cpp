#include <Atoms.h>
#include <Eigen/Dense>
#include <gtest/gtest.h>
#include <iostream>

// Demonstrate some basic assertions.
TEST(VerletTest, numAtoms) {
    Eigen::ArrayXd test_vec = Eigen::ArrayXd((1));
    Eigen::Array3Xd test_mat(3, 1);

    test_mat << 0, 0, 0;
    Atoms atom = Atoms(test_mat, test_vec);
    EXPECT_EQ(atom.nb_atoms(), 1);

    // double force = 1e-9
}