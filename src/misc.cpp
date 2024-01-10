#include "../include/misc.hpp"

#include <fstream>

double u0(double x, double y) {
  if (x < 0) x += 2 * M_PI;  // periodic extension
  if (y < 0) y += 2 * M_PI;  // periodic extension

  return sin(x) + sin(y) + sin(x + y);
}

void write(double *x, int nx, int ny, double t, ofstream &file) {
  file << t << endl;
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny - 1; j++) {
      file << x[i * ny + j] << " ";
    }
    file << x[i * ny + ny - 1] << endl;
  }
  file << endl;
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

void set_wave_numbers(double *kx, double *ky, int nx, int ny_complex) {
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny_complex; j++) {
      if (i <= nx / 2)
        kx[i * ny_complex + j] = i * 1.;
      else
        kx[i * ny_complex + j] = (int)(i - nx) * 1.;
      ky[i * ny_complex + j] = j * 1.;
    }
  }
}

void set_C(double *C, double *kx, double *ky, int nx, int ny_complex, double dt, double nu1, double nu2) {
  double aux;
  for (int i = 0; i < nx * ny_complex; i++) {
    aux = kx[i] * kx[i] + (nu2 / nu1) * ky[i] * ky[i];
    C[i] = 1.0 / (1.0 - dt * aux * (1.0 - nu1 * aux));
  }
}

double lagrange_root(double x2, double dt, double y0, double y1, double y2) {
  // The Lagrange polynomial of degree 2 that passes through the points (x0, y0), (x0 + dt, y1), (x0 + 2 * dt, y2) is:
  // p(x) = 1 / (2 * dt ^ 2) * [ x ^ 2 * (y0 - 2 * y1 + y2) - x * (y0 * (x1 + x2) - 2 * y1 * (x0 + x2) + y2 * (x0 + x1)) + y0 * x1 * x2 - 2 * y1 * x0 * x2 + y2 * x0 * x1 ]
  double x0 = x2 - 2 * dt;
  double x1 = x2 - dt;
  double a = y0 - 2 * y1 + y2;
  double b = -(y0 * (x1 + x2) - 2 * y1 * (x0 + x2) + y2 * (x0 + x1));
  double c = y0 * x1 * x2 - 2 * y1 * x0 * x2 + y2 * x0 * x1;
  // delta is positive if we have entered in this function
  double delta_sq = sqrt(b * b - 4 * a * c);
  double root_2 = (-b - delta_sq) / (2 * a);
  double root_1 = (-b + delta_sq) / (2 * a);
  const double TOL = 1e-10;

  // we check if x0 <= root_1 <= x1 or x1 <= root_1 <= x2
  if ((x0 - TOL < root_1 && root_1 < x1 + TOL) || (x1 - TOL < root_1 && root_1 < x2 + TOL))
    return root_1;
  else
    return root_2;
}

double lagrange_eval(double t, double x2, double dt, double y0, double y1, double y2) {
  double x0 = x2 - 2 * dt;
  double x1 = x2 - dt;
  double a = y0 - 2 * y1 + y2;
  double b = -(y0 * (x1 + x2) - 2 * y1 * (x0 + x2) + y2 * (x0 + x1));
  double c = y0 * x1 * x2 - 2 * y1 * x0 * x2 + y2 * x0 * x1;
  // evaluate the Lagrange polynomial at t
  return (c + t * (b + t * a)) / (2 * dt * dt);
}