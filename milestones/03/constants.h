#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <map>

const double BOLTZMANN = 8.617333262e-5;
// https://github.com/lumol-org/lumol/blob/813ee362428bd19e2797d917d2f32206a83f65d6/lumol-core/src/sys/config/mass.rs
const std::map<std::string, double> ATOMIC_MASSES = {
    {"H", 1.008},        {"He", 4.002602},    {"Li", 6.94},
    {"Be", 9.012182},    {"B", 10.81},        {"C", 12.011},
    {"N", 14.007},       {"O", 15.999},       {"F", 18.9984032},
    {"Ne", 20.1797},     {"Na", 22.98976928}, {"Mg", 24.305},
    {"Al", 26.9815386},  {"Si", 28.085},      {"P", 30.973762},
    {"S", 32.06},        {"Cl", 35.45},       {"Ar", 39.948},
    {"K", 39.0983},      {"Ca", 40.078},      {"Sc", 44.955912},
    {"Ti", 47.867},      {"V", 50.9415},      {"Cr", 51.9961},
    {"Mn", 54.938045},   {"Fe", 55.845},      {"Co", 58.933195},
    {"Ni", 58.6934},     {"Cu", 63.546},      {"Zn", 65.38},
    {"Ga", 69.723},      {"Ge", 72.63},       {"As", 74.9216},
    {"Se", 78.96},       {"Br", 79.904},      {"Kr", 83.798},
    {"Rb", 85.4678},     {"Sr", 87.62},       {"Y", 88.90585},
    {"Zr", 91.224},      {"Nb", 92.90638},    {"Mo", 95.96},
    {"Tc", 97.0},        {"Ru", 101.07},      {"Rh", 102.9055},
    {"Pd", 106.42},      {"Ag", 107.8682},    {"Cd", 112.411},
    {"In", 114.818},     {"Sn", 118.71},      {"Sb", 121.76},
    {"Te", 127.6},       {"I", 126.90447},    {"Xe", 131.293},
    {"Cs", 132.9054519}, {"Ba", 137.327},     {"La", 138.90547},
    {"Ce", 140.116},     {"Pr", 140.90765},   {"Nd", 144.242},
    {"Pm", 145.0},       {"Sm", 150.36},      {"Eu", 151.964},
    {"Gd", 157.25},      {"Tb", 158.92535},   {"Dy", 162.5},
    {"Ho", 164.93032},   {"Er", 167.259},     {"Tm", 168.93421},
    {"Yb", 173.054},     {"Lu", 174.9668},    {"Hf", 178.49},
    {"Ta", 180.94788},   {"W", 183.84},       {"Re", 186.207},
    {"Os", 190.23},      {"Ir", 192.217},     {"Pt", 195.084},
    {"Au", 196.966569},  {"Hg", 200.592},     {"Tl", 204.38},
    {"Pb", 207.2},       {"Bi", 208.9804},    {"Po", 209.0},
    {"At", 210.0},       {"Rn", 222.0},       {"Fr", 223.0},
    {"Ra", 226.0},       {"Ac", 227.0},       {"Th", 232.03806},
    {"Pa", 231.03588},   {"U", 238.02891},    {"Np", 237.0},
    {"Pu", 244.0},       {"Am", 243.0},       {"Cm", 247.0},
    {"Bk", 247.0},       {"Cf", 251.0},       {"Es", 252.0},
    {"Fm", 257.0},       {"Md", 258.0},       {"No", 259.0},
    {"Lr", 262.0},       {"Rf", 267.0},       {"Db", 270.0},
    {"Sg", 271.0},       {"Bh", 270.0},       {"Hs", 277.0},
    {"Mt", 276.0},       {"Ds", 281.0},       {"Rg", 282.0},
    {"Cn", 285.0},       {"Uut", 285.0},      {"Fl", 289.0},
    {"Mc", 289.0},       {"Lv", 293.0},       {"Ts", 294.0},
    {"Og", 294.0},
};
#endif // CONSTANTS_H
