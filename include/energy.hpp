#ifndef ENERGY_HPP
#define ENERGY_HPP

#include <fftw3.h>
typedef struct My_double {
  int int_part;
  double frac_part;
} My_double;
// @brief: computes the L²-energy (squared) of the system, which by definition is the square of the L²-norm of the solution
// @param x: the data of the system at a given time (in the real space).
// @param nx: number of points in the x-direction
// @param ny: number of points in the y-direction
// @return: the L²-energy of the system
double E(double *x, int nx, int ny);
My_double My_E(double *x, int nx, int ny);
double l2(double *x, int nx, int ny);
double E_cplx(fftw_complex *x, int nx, int ny);
#endif  // ENERGY_HPP