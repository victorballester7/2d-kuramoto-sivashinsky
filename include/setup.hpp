#ifndef SETUP_HPP
#define SETUP_HPP

#include <iostream>

#include "LinearAlgebra.hpp"

using namespace std;

// typedef struct args {
//   double nu1, nu2;
//   Vector<int> nn;
//   Vector<double> k1, k2;
//   Vector<double> tmp;  // tmp = (k₁² +(ν₂ / ν₁) * k₂²)
// } Args;

class Args {
 private:
  Vector<double> aux;

 public:
  int nx, ny;
  double h;
  Vector<int> nn;
  Vector<double> k1, k2;  // wave numbers for the frequency space in y as a meshgrid
  // -----------------------------------------------------
  // k1                   k2
  // x0 x0 x0 x0 ... x0   y0 y0 y1 y1 ... ym ym
  // x1 x1 x1 x1 ... x1   y0 y0 y1 y1 ... ym ym
  // .  .  .  .  ... .    .  .  .  .  ... .  .
  // .  .  .  .  ... .    .  .  .  .  ... .  .
  // .  .  .  .  ... .    .  .  .  .  ... .  .
  // xn xn xn xn ... xn   y0 y0 y1 y1 ... ym ym
  double nu1, nu2;
  Vector<double> tmp;  // temporary vector for iterative computations

  Args(int nx_, int ny_, double h_, double nu1_, double nu2_) : nx(nx_), ny(ny_), h(h_), nu1(nu1_), nu2(nu2_), nn(nx_, ny_), k1(2 * nx_ * ny_), k2(2 * nx_ * ny_), tmp(2 * nx_ * ny_), aux(2 * nx_ * ny_) {
    set_wave_numbers();
    aux = k1 * k1 + (nu2 / nu1) * k2 * k2;
    set_tmp();
  }

  void set_wave_numbers() {
    for (int i = 0; i < nx; i++) {
      for (int j = 0; j < 2 * ny; j++) {
        if (i <= nx / 2)
          k1[2 * i * ny + j] = i * 1.;
        else
          k1[2 * i * ny + j] = (i - nx) * 1.;
        if (j <= ny)
          k2[2 * i * ny + j] = (j / 2) * 1.;
        else
          k2[2 * i * ny + j] = (j / 2 - ny) * 1.;
      }
    }
  }

  void set_tmp() {
    tmp = 1.0 / (1.0 + h * aux * ((nu1 * aux) - 1.0));
  }
};

double u0(double x, double y);
double write(Vector<double> &x, Vector<int> nn, double t, ofstream &file, bool print_complex = false);
void stepFiniteDiff(Vector<double> &x, double &t, double h, Args &prm);
void set_data(Vector<double> &x, Args &prm);

#endif  // SETUP_HPP