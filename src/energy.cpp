#include "../include/energy.hpp"

double l2_squared(vector<double> &x, const Args &prm) {
  double E = 0.0;
  for (uint i = 0; i < prm.nx; i++) {
    for (uint j = 0; j < prm.ny; j++) {
      E += x[2 * (i * prm.ny + j)] * x[2 * (i * prm.ny + j)] + x[2 * (i * prm.ny + j) + 1] * x[2 * (i * prm.ny + j) + 1];
    }
  }
  return E;
}

double E(vector<double> &x, const Args &prm) {
  double energy = 0.0;
  for (uint i = 0; i < prm.nx; i++) {
    for (uint j = 0; j < prm.ny; j++) {
      energy += x[2 * (i * prm.ny + j)] * x[2 * (i * prm.ny + j)];
    }
  }
  return energy * (2 * M_PI / prm.nx) * (2 * M_PI / prm.ny);
}

// double dE(vector<double> &x, const Args &prm) {
//   double dE = 0.0;
//   for (uint i = 0; i < prm.nx; i++) {
//     for (uint j = 0; j < prm.ny; j++) {
//       dE += x[2 * (i * prm.ny + j)] * x[2 * (i * prm.ny + j)] + x[2 * (i * prm.ny + j) + 1] * x[2 * (i * prm.ny + j) + 1];
//     }
//   }
//   return dE;
// }