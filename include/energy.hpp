#ifndef ENERGY_HPP
#define ENERGY_HPP

// @brief: computes the L²-energy (squared) of the system, which by definition is the square of the L²-norm of the solution
// @param x: the data of the system at a given time (in the real space).
// @param nx: number of points in the x-direction
// @param ny: number of points in the y-direction
// @return: the L²-energy of the system
double E(double *x, int nx, int ny);

#endif  // ENERGY_HPP