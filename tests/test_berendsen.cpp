#include <gtest/gtest.h>

#include "Atoms.h"
#include "ForcesEnergies.h"
#include "LatticeGen.h"
#include "berendsen.h"
#include "verlet.h"

TEST(BerendsenTest, CorrectRescaling) {
    // construct fresh system, T0 ~ 0.55
    Atoms system = cubic(50, 1.0, 1.008, 0.1);
    // simu parameters
    double timestep = 0.01;
    int steps = 50;
    // berendsen parameters
    double goal_temp = 0.2;
    double tau = 1;
    // run some steps to make sure velocities are there and
    // the "true" starting temperature is reached
    for (int i = 0; i < steps; i++) {
        verlet1(system, timestep);
        lj_direct_summation(system, 1.0, 1.0);
        verlet2(system, timestep);
    }
    // test starts here
    double t_before = temperature_lj(system);
    // apply thermostat
    berendsen_thermostat_lj(system, goal_temp, timestep, tau);
    double t_after = temperature_lj(system);
    // expected value
    double t_expected =
        goal_temp + (t_before - goal_temp) * exp(-(timestep / tau));

    // check if calculated solution is same as simulated
    EXPECT_NEAR(t_after, t_expected, 1e-5);
}
