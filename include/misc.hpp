#ifndef SETUP_HPP
#define SETUP_HPP

#include <iostream>

#include "vector.hpp"

using namespace std;

// typedef struct args {
//   double nu1, nu2;
//   vector<uint> nn;
//   vector<double> k1, k2;
//   vector<double> tmp;  // tmp = (k₁² +(ν₂ / ν₁) * k₂²)
// } Args;

class Args {
 public:
  uint nx, ny;
  double h;
  vector<uint> nn = {nx, ny};  // number of points in each direction
  vector<double> k1, k2;       // wave numbers for the frequency space in y as a meshgrid
  // -----------------------------------------------------
  // k1                   k2
  // x0 x0 x0 x0 ... x0   y0 y0 y1 y1 ... ym ym
  // x1 x1 x1 x1 ... x1   y0 y0 y1 y1 ... ym ym
  // .  .  .  .  ... .    .  .  .  .  ... .  .
  // .  .  .  .  ... .    .  .  .  .  ... .  .
  // .  .  .  .  ... .    .  .  .  .  ... .  .
  // xn xn xn xn ... xn   y0 y0 y1 y1 ... ym ym
  double nu1, nu2;
  vector<double> tmp;  // temporary vector for iterative computations
  vector<double> L;    // linear part of the equation û_t = L * û

  Args(uint nx_, uint ny_, double h_, double nu1_, double nu2_) : nx(nx_), ny(ny_), h(h_), k1(2 * nx_ * ny_), k2(2 * nx_ * ny_), nu1(nu1_), nu2(nu2_), tmp(2 * nx_ * ny_), L(2 * nx_ * ny_), aux(2 * nx_ * ny_) {
    set_wave_numbers();
    aux = (k1 * k1) + ((nu2 / nu1) * k2 * k2);
    L = aux * (1.0 - (nu1 * aux));
    set_tmp();
  }

  void set_wave_numbers() {
    for (uint i = 0; i < nx; i++) {
      for (uint j = 0; j < 2 * ny; j++) {
        if (i <= nx / 2)
          k1[2 * i * ny + j] = i * 1.;
        else
          k1[2 * i * ny + j] = (int)(i - nx) * 1.;
        if (j <= ny + 1)
          k2[2 * i * ny + j] = (int)(j / 2) * 1.;
        else
          k2[2 * i * ny + j] = (int)(j / 2 - ny) * 1.;
      }
    }
  }

  void set_tmp() {
    tmp = 1.0 / (1.0 - (h * L));
  }

 private:
  vector<double> aux;
};

// @brief: computes the initial condition
// @param x: x-coordinate
// @param y: y-coordinate
// @return: the value of the initial condition at (x, y)
double u0(double x, double y);

// @brief: writes the data to a file
// @param x: the data to be written
// @param nn: number of points in each direction
// @param t: time at which the data is written
// @param file: file to write the data to
// @param print_complex: whether to print the imaginary part of the data or not
void write(vector<double> &x, vector<uint> nn, double t, ofstream &file);
void write_old(vector<double> &x, vector<uint> nn, double t, ofstream &file, bool print_complex = false);

// @brief: computes the mean of the data as the double integral of the data over the domain divided by the area of the domain
// @param x: the data
// @param prm: parameters of the system
// @return: the mean
double mean(const vector<double> &x, const Args &prm);

// @brief: computes the initial condition on the whole vector
// @param x: the vector to be updated
// @param prm: parameters of the system
void set_data(vector<double> &x, const Args &prm);

#endif  // SETUP_HPP