#ifndef ENERGY_HPP
#define ENERGY_HPP

// @brief: computes the L²-energy of the system, which by definition is the square of the L²-norm of the solution
// @param x: the data of the system at a given time (in the real space). Note that it contains both the real and imaginary part, but we only use the real part.
// @param prm: parameters of the system
// @return: the L²-energy of the system
// double l2_squared(vector<double> &x, const Args &prm);
double E(double *x, int nx, int ny);

#endif  // ENERGY_HPP