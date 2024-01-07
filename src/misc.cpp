#include "../include/misc.hpp"

#include <fstream>

#include "../include/fft.hpp"
// #include "../include/matplotlibcpp.h"
#include "../include/vector.hpp"

// namespace plt = matplotlibcpp;

double u0(double x, double y) {  // initial condition
  if (x < 0) x += 2 * M_PI;      // periodic extension
  if (y < 0) y += 2 * M_PI;      // periodic extension

  // return sin(x) + sin(2 * y) + cos(7 * exp(-x));
  return sin(x) + sin(y) + sin(x + y);
}

// double heat_kernel(double t, double x, double y) {
//   return 1.0 / (4 * M_PI * t) * exp(-(x * x + y * y) / (4 * t));
// }

// double u(double t, double x, double y, int nx, int ny) {
//   // integrate numerically f * H, where H is the heat kernel
//   double xmin = 0, xmax = 2 * M_PI;
//   double ymin = 0, ymax = 2 * M_PI;
//   double hx = (xmax - xmin) / nx;
//   double hy = (ymax - ymin) / ny;
//   double sum = 0.0;

//   for (double i = xmin; i < xmax; i += hx) {
//     for (double j = ymin; j < ymax; j += hy) {
//       sum += f(i, j) * heat_kernel(t, x - i, y - j);
//     }
//   }
//   return sum * hx * hy;
// }

void write(vector<double> &x, vector<uint> nn, double t, ofstream &file) {  // write into file
  file << t << endl;
  for (uint i = 0; i < nn[0]; i++) {
    for (uint j = 0; j < nn[1] - 1; j++) {
      file << x[2 * (i * nn[1] + j)] << " ";
    }
    file << x[2 * (i * nn[1] + nn[1] - 1)] << endl;
  }
  file << endl;
}

void write_old(vector<double> &x, vector<uint> nn, double t, ofstream &file, bool print_complex) {  // write into file
  file << t << endl;
#define round_to_zero(x) (fabs(x) < 1e-10 ? 0 : x)

  for (uint i = 0; i < nn[0]; i++) {
    for (uint j = 0; j < 2 * nn[1] - 2; j += 2) {
      if (print_complex)
        file << round_to_zero(x[2 * i * nn[1] + j]) << showpos << round_to_zero(x[2 * i * nn[1] + j + 1]) << noshowpos << "i ";
      else
        file << round_to_zero(x[2 * i * nn[1] + j]) << " ";
    }
    if (print_complex)
      file << round_to_zero(x[2 * i * nn[1] + 2 * nn[1] - 2]) << showpos << round_to_zero(x[2 * i * nn[1] + 2 * nn[1] - 1]) << noshowpos << "i" << endl;
    else
      file << round_to_zero(x[2 * i * nn[1] + 2 * nn[1] - 2]) << endl;
  }
  file << endl;
#undef round_to_zero
}

// double mean(const vector<double> &x, const Args &prm) {
//   double sum = 0.0;
//   for (uint i = 0; i < prm.nx; i++) {
//     for (uint j = 0; j < 2 * prm.ny; j += 2) {
//       sum += x[2 * i * prm.ny + j];
//     }
//   }
//   return sum / (prm.nx * prm.ny);  // = 1 / (2 * M_PI) ^ 2 * sum * (2 * M_PI / prm.nx) * (2 * M_PI / prm.ny)
// }

void set_data(vector<double> &x, const Args &prm) {
  double a1 = 2 * M_PI / prm.nx;
  double a2 = 2 * M_PI / prm.ny;
  for (uint i = 0; i < prm.nx; i++) {
    for (uint j = 0; j < 2 * prm.ny; j += 2) {
      x[2 * i * prm.ny + j] = u0(prm.k1[2 * i * prm.ny + j] * a1, prm.k2[2 * i * prm.ny + j] * a2);
      x[2 * i * prm.ny + j + 1] = 0;
    }
  }
}
