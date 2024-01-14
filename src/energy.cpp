#include "../include/energy.hpp"

#include <cmath>

double E(double *x, int nx, int ny) {
  double energy = 0.0;
  for (int i = 0; i < nx * ny; i++) {
    energy += x[i] * x[i];
  }
  return energy * (2 * M_PI / nx) * (2 * M_PI / ny);
}

My_double My_E(double *x, int nx, int ny) {
  My_double energy = {0, 0.0};
  // energy.int_part will be the sum of the integer part and the first 3 digits of the fractional part of athe square of the terms in the array
  // energy.frac_part will be the sum of the remaining digits of the fractional part of the square of the terms in the array
  int int_part = 0;
  for (int i = 0; i < nx * ny; i++) {
    int_part = (int)(x[i] * x[i]);
    energy.frac_part += x[i] * x[i] - int_part;
    energy.int_part += int_part;
    // int_part = x[i];
    // frac_part = x[i] - int_part;
    // energy.int_part += int_part * int_part;
    // energy.frac_part += frac_part * (frac_part + 2 * int_part);
    // tmp = energy.frac_part;
    // energy.int_part += tmp;
    // energy.frac_part -= tmp;
    // energy.frac_part += x[i] * x[i];
  }
  // energy.int_part *= (2 * M_PI / nx) * (2 * M_PI / ny);
  // energy.frac_part *= (2 * M_PI / nx) * (2 * M_PI / ny);
  return energy;
}

double l2(double *x, int nx, int ny) {
  double l = 0.0;
  for (int i = 0; i < nx * ny; i++) {
    l += x[i] * x[i];
  }
  return sqrt(l / (nx * ny));
}

double E_cplx(fftw_complex *x, int nx, int ny) {
  double energy = 0.0;
  for (int i = 0; i < nx * ny; i++) {
    energy += (x[i][0] * x[i][0] + x[i][1] * x[i][1]);
  }
  return energy * (2 * M_PI / nx) * (2 * M_PI / ny);
}