#include "../include/energy.hpp"

// double l2_squared(vector<double> &x, const Args &prm) {
//   double E = 0.0;
//   for (int i = 0; i < prm.nx; i++) {
//     for (int j = 0; j < prm.ny; j++) {
//       E += x[2 * (i * prm.ny + j)] * x[2 * (i * prm.ny + j)] + x[2 * (i * prm.ny + j) + 1] * x[2 * (i * prm.ny + j) + 1];
//     }
//   }
//   return E;
// }

double E(double *x, int nx, int ny) {
  double energy = 0.0;
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      energy += x[i * ny + j] * x[i * ny + j];
    }
  }
  return energy * (2 * M_PI / nx) * (2 * M_PI / ny);
}

// double dE(vector<double> &x, const Args &prm) {
//   double dE = 0.0;
//   for (int i = 0; i < prm.nx; i++) {
//     for (int j = 0; j < prm.ny; j++) {
//       dE += x[2 * (i * prm.ny + j)] * x[2 * (i * prm.ny + j)] + x[2 * (i * prm.ny + j) + 1] * x[2 * (i * prm.ny + j) + 1];
//     }
//   }
//   return dE;
// }