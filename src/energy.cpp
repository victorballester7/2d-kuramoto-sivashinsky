#include "../include/energy.hpp"

#include <cmath>

double E(double *x, int nx, int ny) {
  double energy = 0.0;
  for (int i = 0; i < nx * ny; i++) {
    energy += x[i] * x[i];
  }
  return energy * (2 * M_PI / nx) * (2 * M_PI / ny);
}

double l2(double *x, int nx, int ny) {
  double l = 0.0;
  for (int i = 0; i < nx * ny; i++) {
    l += x[i] * x[i];
  }
  return sqrt(l);
}

double E_cplx(fftw_complex *x, int nx, int ny) {
  double energy = 0.0;
  for (int i = 0; i < nx * ny; i++) {
    energy += (x[i][0] * x[i][0] + x[i][1] * x[i][1]);
  }
  return energy * (2 * M_PI / nx) * (2 * M_PI / ny);
}