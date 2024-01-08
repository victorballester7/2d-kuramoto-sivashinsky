#include <cmath>

double E(double *x, int nx, int ny) {
  double energy = 0.0;
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      energy += x[i * ny + j] * x[i * ny + j];
    }
  }
  return energy * (2 * M_PI / nx) * (2 * M_PI / ny);
}