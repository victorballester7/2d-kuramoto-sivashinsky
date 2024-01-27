#include "../include/misc.hpp"

#include <fftw3.h>
#include <math.h>

#include <fstream>

double u0(double x, double y) {
  if (x < 0) x += 2 * M_PI;  // periodic extension
  if (y < 0) y += 2 * M_PI;  // periodic extension

  return sin(x) + sin(y) + sin(x + y);
  // return sin(x) + sin(y) + cos(x + y) + sin(4 * x + 4 * y) + cos(7 * x) + cos(7 * y);
}

double U(double x, double y, double t) {
  if (x < 0) x += 2 * M_PI;  // periodic extension
  if (y < 0) y += 2 * M_PI;  // periodic extension

  return sin(x) + sin(y) + sin(x + y) + sin(t);
}

double g(double x, double y, double t) {
  if (x < 0) x += 2 * M_PI;  // periodic extension
  if (y < 0) y += 2 * M_PI;  // periodic extension

  // return 1 + cos(2 * x + 2 * y) / 2 + cos(2 * x) / 4 + cos(2 * y) / 4 + cos(y + 2 * x) / 2 + cos(y) / 2 + cos(x + 2 * y) / 2 + cos(x) / 2 + 2 * sin(x + y) + cos(t);
  return 1 + cos(x + 2 * y) / 2 + cos(y + 2 * x) / 2 + cos(2 * x + 2 * y) / 2 + cos(2 * y) / 4 + cos(2 * x) / 4 + cos(x) / 2 + cos(y) / 2 + cos(t) + (6 * sin(x + y)) / 5 - sin(x) / 5 - sin(y) / 5;
  // return cos(t) + (6 * sin(x + y)) / 5 - sin(x) / 5 - sin(y) / 5;
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

void write_fourier(fftw_complex *x, int nx, int ny_complex, double t, ofstream &file) {
  file << t << endl;
  for (int i = 0; i < nx / 2 + 1; i++) {
    for (int j = 0; j < ny_complex - 1; j++) {
      if (i == 0 && j == 0) {
        file << 0 << " ";  // we set the mean to 0 (without editing the data in the vector)
        continue;
      }
      file << x[i * ny_complex + j][0] * x[i * ny_complex + j][0] + x[i * ny_complex + j][1] * x[i * ny_complex + j][1] << " ";
    }
    file << x[i * ny_complex + ny_complex - 1][0] * x[i * ny_complex + ny_complex - 1][0] + x[i * ny_complex + ny_complex - 1][1] * x[i * ny_complex + ny_complex - 1][1] << endl;
  }
  file << endl;
}

void write_bis(fftw_complex *x, int nx, int ny_complex, double t, ofstream &file) {
#define round_to_zero(y) (fabs(y) < 1e-10 ? 0 : y)
  file << t << endl;
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny_complex - 1; j++) {
      file << round_to_zero(x[i * ny_complex + j][0]) << " " << round_to_zero(x[i * ny_complex + j][1]) << " | ";
    }
    file << round_to_zero(x[i * ny_complex + ny_complex - 1][0]) << " " << round_to_zero(x[i * ny_complex + ny_complex - 1][1]) << endl;
  }
  file << endl;
#undef round_to_zero
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

void set_C1(double *C1, double *kx, double *ky, int nx, int ny_complex, double dt, double nu1, double nu2, double c) {
  double aux;
  for (int i = 0; i < nx * ny_complex; i++) {
    aux = kx[i] * kx[i] + (nu2 / nu1) * ky[i] * ky[i];
    C1[i] = 1.0 / (1.0 + dt * (aux * (nu1 * aux - 1.0) + c));
  }
}

void set_C2(double *C2, double *kx, double *ky, int nx, int ny_complex, double dt, double nu1, double nu2, double c) {
  double aux;
  for (int i = 0; i < nx * ny_complex; i++) {
    aux = kx[i] * kx[i] + (nu2 / nu1) * ky[i] * ky[i];
    C2[i] = 1.0 / (1.5 + dt * (aux * (nu1 * aux - 1.0) + c));
  }
}

double max(double *x, uint dim) {
  double max = fabs(x[0]);
  double tmp;
  for (uint i = 1; i < dim; i++) {
    tmp = fabs(x[i]);
    if (tmp > max) max = tmp;
  }
  return max;
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

double lagrange_root(double x2, double dt, My_double y0, My_double y1, My_double y2) {
  // The Lagrange polynomial of degree 2 that passes through the points (x2 - 2 * dt, y0), (x2 - dt, y1), (x2, y2) is:
  // p(x) = 1 / (2 * dt ^ 2) * [ x ^ 2 * (y0 - 2 * y1 + y2) - x * (y0 * (x1 + x2) - 2 * y1 * (x0 + x2) + y2 * (x0 + x1)) + y0 * x1 * x2 - 2 * y1 * x0 * x2 + y2 * x0 * x1 ]
  double a = (y0.int_part - 2 * y1.int_part + y2.int_part) + (y0.frac_part - 2 * y1.frac_part + y2.frac_part);
  double tmp = (y0.int_part - 4 * y1.int_part + 3 * y2.int_part) + (y0.frac_part - 4 * y1.frac_part + 3 * y2.frac_part);
  double b = -2 * a * x2 + tmp * dt;
  double c = x2 * (a * x2 - tmp * dt) + 2 * (y2.int_part + y2.frac_part) * dt * dt;
  // delta is positive if we have entered in this function
  double delta_sq = sqrt(b * b - 4 * a * c);

  double root_2 = (b > 0) ? (-b - delta_sq) / (2 * a) : (2 * c) / (-b + delta_sq);  // in order to avoid numerical round-off errors
  double root_1 = (b < 0) ? (-b + delta_sq) / (2 * a) : (2 * c) / (-b - delta_sq);  // in order to avoid numerical round-off errors
  const double TOL = 1e-10;

  // we check if x0 <= root_1 <= x1 or x1 <= root_1 <= x2
  if ((x2 - 2 * dt - TOL < root_1 && root_1 < x2 - dt + TOL) || (x2 - dt - TOL < root_1 && root_1 < x2 + TOL))
    return root_1;
  else
    return root_2;
}

double lagrange_eval(double t, double x2, double dt, My_double y0, My_double y1, My_double y2) {
  double a = (y0.int_part - 2 * y1.int_part + y2.int_part) + (y0.frac_part - 2 * y1.frac_part + y2.frac_part);
  double tmp = (y0.int_part - 4 * y1.int_part + 3 * y2.int_part) + (y0.frac_part - 4 * y1.frac_part + 3 * y2.frac_part);
  double b = -2 * a * x2 + tmp * dt;
  double c = x2 * (a * x2 - tmp * dt) + 2 * (y2.int_part + y2.frac_part) * dt * dt;
  // evaluate the Lagrange polynomial at t
  return (c + t * (b + t * a)) / (2 * dt * dt);
}

void print(string str, int64_t time) {
  // default time is in microseconds
  string microsec = "\u03BCs";
  string millisec = "ms";
  string sec = "s";
  string min = "min";
  string hour = "h";
  double time_d;
  if (time > 10000) {
    time /= 1000;  // now time is in milliseconds
    if (time > 10000) {
      time /= 1000;  // now time is in seconds
      if (time > 100) {
        time_d = (double)time;
        time_d /= 60;  // now time is in minutes
        if (time_d > 100) {
          time_d /= 60;  // now time is in hours
          cout << str << time_d << " " << hour << endl;
          return;
        }
        cout << str << time_d << " " << min << endl;
        return;
      }
      cout << str << time << " " << sec << endl;
      return;
    }
    cout << str << time << " " << millisec << endl;
    return;
  }
  cout << str << time << " " << microsec << endl;
  return;
}