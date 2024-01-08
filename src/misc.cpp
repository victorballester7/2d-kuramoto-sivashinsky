#include "../include/misc.hpp"

#include <fstream>

// #include "../include/fft.hpp"
// #include "../include/matplotlibcpp.h"
#include "../include/vector.hpp"

double u0(double x, double y) {  // initial condition
  if (x < 0) x += 2 * M_PI;      // periodic extension
  if (y < 0) y += 2 * M_PI;      // periodic extension

  return sin(x) + sin(y) + sin(x + y);
}

// double d_y_u0(double x, double y) {  // derivative of the initial condition with respect to y
//   if (x < 0) x += 2 * M_PI;          // periodic extension
//   if (y < 0) y += 2 * M_PI;          // periodic extension

//   return cos(y) + cos(x + y) + 4 * sin(4 * x) * cos(4 * y) + 7 * cos(5 * x) * cos(7 * y) + 8. * cos(3 * (x + y));
// }

// double d_x_u0(double x, double y) {  // derivative of the initial condition with respect to x
//   if (x < 0) x += 2 * M_PI;          // periodic extension
//   if (y < 0) y += 2 * M_PI;          // periodic extension

//   return cos(x) + cos(x + y) + 4 * cos(4 * x) * sin(4 * y) - 7 * sin(7 * x) - 5 * sin(5 * x) * sin(9 * y) + 2 * cos(2 * x) + 8. * cos(3 * (x + y)) + 7. * cos(7 * x) + 2 * sin(6 * x);
// }

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

void write(double *x, int nx, int ny, double t, ofstream &file) {  // write into file
  file << t << endl;
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny - 1; j++) {
      file << x[i * ny + j] << " ";
    }
    file << x[i * ny + ny - 1] << endl;
  }
  file << endl;
}

// void write_old(fftw_complex *x, int nx, int ny, double t, ofstream &file) {  // write into file
//   file << t << endl;
// #define round_to_zero(y) (fabs(y) < 1e-10 ? 0 : y)

//   for (int i = 0; i < nx; i++) {
//     for (int j = 0; j < ny - 1; j++) {
//       file << round_to_zero(x[i * ny + j][0]) << " " << round_to_zero(x[i * ny + j][1]) << noshowpos << " | ";
//     }
//     file << round_to_zero(x[i * ny + ny - 1][0]) << " " << round_to_zero(x[i * ny + ny - 1][1]) << noshowpos << " | " << endl;
//   }
//   file << endl;
// #undef round_to_zero
// }

// void write_old_old(vector<double> &x, vector<uint> nn, double t, ofstream &file, bool print_complex) {  // write into file
//   file << t << endl;
// #define round_to_zero(x) (fabs(x) < 1e-10 ? 0 : x)

//   for (int i = 0; i < nn[0]; i++) {
//     for (int j = 0; j < 2 * nn[1] - 2; j += 2) {
//       if (print_complex)
//         file << round_to_zero(x[2 * i * nn[1] + j]) << showpos << round_to_zero(x[2 * i * nn[1] + j + 1]) << noshowpos << "i ";
//       else
//         file << round_to_zero(x[2 * i * nn[1] + j]) << " ";
//     }
//     if (print_complex)
//       file << round_to_zero(x[2 * i * nn[1] + 2 * nn[1] - 2]) << showpos << round_to_zero(x[2 * i * nn[1] + 2 * nn[1] - 1]) << noshowpos << "i" << endl;
//     else
//       file << round_to_zero(x[2 * i * nn[1] + 2 * nn[1] - 2]) << endl;
//   }
//   file << endl;
// #undef round_to_zero
// }

// double mean(const vector<double> &x, const Args &prm) {
//   double sum = 0.0;
//   for (int i = 0; i < prm.nx; i++) {
//     for (int j = 0; j < 2 * prm.ny; j += 2) {
//       sum += x[2 * i * prm.ny + j];
//     }
//   }
//   return sum / (prm.nx * prm.ny);  // = 1 / (2 * M_PI) ^ 2 * sum * (2 * M_PI / prm.nx) * (2 * M_PI / prm.ny)
// }

double mean(double *x, int nx, int ny) {
  double sum = 0.0;
  for (int i = 0; i < nx * ny; i++) {
    sum += x[i];
  }
  return sum / (nx * ny);  // = 1 / (2 * M_PI) ^ 2 * sum * (2 * M_PI / prm.nx) * (2 * M_PI / prm.ny)
}

void set_data(double *x, int nx, int ny) {
  double a1 = 2 * M_PI / nx;
  double a2 = 2 * M_PI / ny;
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      x[i * ny + j] = u0(i * a1, j * a2);
    }
  }
}
void set_wave_numbers(double *k1, double *k2, int nx, int ny_complex) {
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny_complex; j++) {
      if (i <= nx / 2)
        k1[i * ny_complex + j] = i * 1.;
      else
        k1[i * ny_complex + j] = (int)(i - nx) * 1.;
      k2[i * ny_complex + j] = j * 1.;
    }
  }
}

void set_C(double *tmp, double *k1, double *k2, int nx, int ny_complex, double dt, double nu1, double nu2) {
  double aux;
  for (int i = 0; i < nx * ny_complex; i++) {
    aux = k1[i] * k1[i] + (nu2 / nu1) * k2[i] * k2[i];
    tmp[i] = 1.0 / (1.0 - dt * aux * (1.0 - nu1 * aux));
  }
}