#ifndef ENERGY_HPP
#define ENERGY_HPP
#include "misc.hpp"
#include "vector.hpp"

// @brief: computes the L²-energy of the system, which by definition is the square of the L²-norm of the solution
// @param x: the data of the system at a given time (in the real space). Note that it contains both the real and imaginary part, but we only use the real part.
// @param prm: parameters of the system
// @return: the L²-energy of the system
double l2_squared(vector<double> &x, const Args &prm);
double E(vector<double> &x, const Args &prm);
// @brief: computes the derivative of the L²-energy of the system
// @param x: the data of the system at a given time (in the real space). Note that it contains both the real and imaginary part, but we only use the real part.
// @param prm: parameters of the system
// @return: the derivative of the L²-energy of the system
double dE(vector<double> &x, const Args &prm);

#endif  // ENERGY_HPP