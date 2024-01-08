#ifndef MISC_HPP
#define MISC_HPP

#include <complex.h>
#include <fftw3.h>

// #include <complex>
#include <iostream>

// #include "fftw3.h"
#include "vector.hpp"
using namespace std;

// typedef struct args {
//   double nu1, nu2;
//   vector<uint> nn;
//   vector<double> k1, k2;
//   vector<double> tmp;  // tmp = (k₁² +(ν₂ / ν₁) * k₂²)
// } Args;

// class Args {
//  public:
//   int nx, ny;
//   double h;
//   vector<int> nn = {nx, ny};  // number of points in each direction
//   vector<double> k1, k2;       // wave numbers for the frequency space in y as a meshgrid
//   // -----------------------------------------------------
//   // k1                   k2
//   // x0 x0 x0 x0 ... x0   y0 y0 y1 y1 ... ym ym
//   // x1 x1 x1 x1 ... x1   y0 y0 y1 y1 ... ym ym
//   // .  .  .  .  ... .    .  .  .  .  ... .  .
//   // .  .  .  .  ... .    .  .  .  .  ... .  .
//   // .  .  .  .  ... .    .  .  .  .  ... .  .
//   // xn xn xn xn ... xn   y0 y0 y1 y1 ... ym ym
//   double nu1, nu2;
//   vector<double> tmp;  // temporary vector for iterative computations
//   vector<double> L;    // linear part of the equation û_t = L * û

//   Args(uint nx_, uint ny_, double h_, double nu1_, double nu2_) : nx(nx_), ny(ny_), h(h_), k1(nx_ * ny_), k2(nx_ * ny_), nu1(nu1_), nu2(nu2_), tmp(nx_ * ny_), L(nx_ * ny_), aux(nx_ * ny_) {
//     set_wave_numbers();
//     aux = (k1 * k1) + ((nu2 / nu1) * k2 * k2);
//     L = aux * (1.0 - (nu1 * aux));
//     set_tmp();
//   }

//   void set_wave_numbers() {
//     for (int i = 0; i < nx; i++) {
//       for (int j = 0; j < ny; j++) {
//         if (i <= nx / 2)
//           k1[i * ny + j] = i * 1.;
//         else
//           k1[i * ny + j] = (int)(i - nx) * 1.;
//         if (j <= ny + 1)
//           k2[i * ny + j] = j * 1.;
//         else
//           k2[i * ny + j] = (int)(j - ny) * 1.;
//       }
//     }
//   }

//   void set_tmp() {
//     tmp = 1.0 / (1.0 - (h * L));
//   }

//  private:
//   vector<double> aux;
// };

// @brief: computes the initial condition
// @param x: x-coordinate
// @param y: y-coordinate
// @return: the value of the initial condition at (x, y)
double u0(double x, double y);
double d_y_u0(double x, double y);
double d_x_u0(double x, double y);
// @brief: writes the data to a file
// @param x: the data to be written
// @param nn: number of points in each direction
// @param t: time at which the data is written
// @param file: file to write the data to
// @param print_complex: whether to print the imaginary part of the data or not
void write(double *x, int nx, int ny, double t, ofstream &file);

// void write_old(fftw_complex *x, int nx, int ny, double t, ofstream &file);
// void write_old_old(vector<double> &x, vector<uint> nn, double t, ofstream &file, bool print_complex);
// @brief: computes the mean of the data as the double integral of the data over the domain divided by the area of the domain
// @param x: the data
// @param prm: parameters of the system
// @return: the mean
// double mean(const vector<double> &x, const Args &prm);
double mean(double *x, int nx, int ny);
// @brief: computes the initial condition on the whole vector
// @param x: the vector to be updated
// @param prm: parameters of the system
void set_data(double *x, int nx, int ny);

void set_wave_numbers(double *k1, double *k2, int nx, int ny_complex);
void set_C(double *tmp, double *k1, double *k2, int nx, int ny_complex, double dt, double nu1, double nu2);
#endif  // MISC_HPP