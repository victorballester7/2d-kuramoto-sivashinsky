#include <chrono>
#include <fstream>
#include <iostream>

#include "../include/LinearAlgebra.hpp"
#include "../include/fft.hpp"
#include "../include/field.hpp"
#include "../include/rk78.h"

#define UNIT_SECONDS microseconds
#define LABEL_SECONDS "μs"

using namespace std;

double f(double x, double y) {  // initial condition
  // return sin(2 * x) + cos(8 * y);  // should be 2pi-periodic in x and y
  return sin(x);
}

double write(Vector<double> &x, Vector<int> nn, double t, ofstream &file) {  // write into file
  file << t << endl;
  for (int i = 0; i < nn[0]; i++) {
    for (int j = 0; j < 2 * nn[1] - 2; j += 2) {
      file << x[2 * i * nn[1] + j] << " ";
    }
    file << x[2 * i * nn[1] + 2 * nn[1] - 2] << endl;
  }
  file << endl;
  return 0;
}

int main(void) {
  const int nx = 16, ny = 16;                 // number of points in x and y (must be a power of 2)
  const int dim = 2 * nx * ny;                // dimension of the system (the '2' is for the real and imaginary part)
  double t = 0;                               // initial time of integration
  const double T = 5;                         // final time of integration
  double h = 0.1;                             // initial step size
  const double hmin = 0.01;                   // minimum step size
  const double hmax = 0.2;                    // maximum step size
  const double tol = 1e-8;                    // tolerance
  const int maxNumSteps = (int)T / hmin + 1;  // maximum number of steps
  int numSteps = 0;                           // number of steps
  const int maxRepetitionsRK78 = 50;          // maximum number of repetitions of the RK78 algorithm allowed before throwing an error
  const double twoPI = 2 * M_PI;              // 2 * pi
  Vector<double> x(dim);                      // initial condition
  Vector<double> aux(dim);
  Args prm;  // parameters for the field function
  prm.nu1 = 1;
  prm.nu2 = 1;
  prm.k1 = Vector<double>(dim);  // wave numbers for the frequency space in x as a meshgrid
  prm.k2 = Vector<double>(dim);  // wave numbers for the frequency space in y as a meshgrid
  // -----------------------------------------------------
  // k1                   k2
  // x0 x0 x0 x0 ... x0   y0 y0 y1 y1 ... ym ym
  // x1 x1 x1 x1 ... x1   y0 y0 y1 y1 ... ym ym
  // .  .  .  .  ... .    .  .  .  .  ... .  .
  // .  .  .  .  ... .    .  .  .  .  ... .  .
  // .  .  .  .  ... .    .  .  .  .  ... .  .
  // xn xn xn xn ... xn   y0 y0 y1 y1 ... ym ym

  prm.nn = Vector<int>(nx, ny);  // number of points in x and y (must be a power of 2)
  chrono::steady_clock::time_point begin_write, end_write;
  chrono::steady_clock::time_point begin_computations, end_computations;
  int64_t total_write = 0, total_computations = 0;
  begin_computations = chrono::steady_clock::now();
  // set wave numbers
  double a1 = twoPI / nx;
  double a2 = twoPI / ny;
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < 2 * ny; j++) {
      if (i <= nx / 2)
        prm.k1[2 * i * ny + j] = i * a1;
      else
        prm.k1[2 * i * ny + j] = (i - nx) * a1;
      if (j <= ny + 1)
        prm.k2[2 * i * ny + j] = (j / 2) * a2;
      else
        prm.k2[2 * i * ny + j] = (j / 2 - ny) * a2;
    }
  }

  // file handling
  ofstream file;  // output file
  file.open("data/data.txt");

  // set initial condition
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < 2 * ny; j += 2) {
      x[2 * i * ny + j] = f(prm.k1[2 * i * ny + j], prm.k2[2 * i * ny + j]);
      x[2 * i * ny + j + 1] = 0;
    }
  }
  end_computations = chrono::steady_clock::now();
  total_computations += chrono::duration_cast<chrono::UNIT_SECONDS>(end_computations - begin_computations).count();

  // write into file
  begin_write = chrono::steady_clock::now();
  write(x, prm.nn, t, file);
  end_write = chrono::steady_clock::now();
  total_write += chrono::duration_cast<chrono::UNIT_SECONDS>(end_write - begin_write).count();

  // compute the FFT of x (with void fft_n(Vector<double> &data, Vector<int> &nn, const int isign, const int normalize = 1);)

  begin_computations = chrono::steady_clock::now();
  fft_n(x, prm.nn, 1, 1);

  prm.tmp = Vector<double>(dim);  // tmp = k1 ^ 2 + (nu2 / nu1) * k2 ^ 2
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      prm.tmp[i * ny + j] = prm.k1[i] * prm.k1[i] + prm.k2[j] * prm.k2[j] * prm.nu2 / prm.nu1;
    }
  }
  end_computations = chrono::steady_clock::now();
  total_computations += chrono::duration_cast<chrono::UNIT_SECONDS>(end_computations - begin_computations).count();

  do {
    begin_computations = chrono::steady_clock::now();
    if (numSteps == 178) {
      std::cout << "hereeeeeeeeeeeeeeeeeeeeeeeeeeeee: t = " << t << std::endl;
    }
    if (rk78(&t, &x[0], &h, hmin, hmax, tol, maxRepetitionsRK78, dim, field, &prm)) {
      throw std::runtime_error("Error in rk78.");
    }
    // rk78(&t, &x[0], &h, hmin, hmax, tol, dim, field, &prm);
    if (numSteps == 178) {
      std::cout << "hereeeeeeeeeeeeeeeeeeeeeeeeeeeee 2: t = " << t << std::endl;
    }
    // write into file (but only the antitransformed solution)

    aux = x;  // copy x to aux
    fft_n(aux, prm.nn, -1, 1);
    end_computations = chrono::steady_clock::now();
    total_computations += chrono::duration_cast<chrono::UNIT_SECONDS>(end_computations - begin_computations).count();

    begin_write = chrono::steady_clock::now();
    write(aux, prm.nn, t, file);
    end_write = chrono::steady_clock::now();
    total_write += chrono::duration_cast<chrono::UNIT_SECONDS>(end_write - begin_write).count();
    numSteps++;
    // cout << "Step: " << numSteps << endl;
  } while (t < T || numSteps >= maxNumSteps);

  cout << "Total time for computations: " << total_computations << " " << LABEL_SECONDS << endl;
  cout << "Total time for writing: " << total_write << " " << LABEL_SECONDS << endl;

  file.close();
  return 0;
}

// compilate with (the include files are in include/):
// g++ -std=c++11 -I include/ src/main.cpp -o main