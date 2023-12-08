#include <chrono>
#include <fstream>
#include <iostream>

#include "../include/LinearAlgebra.hpp"
#include "../include/fft.hpp"
#include "../include/setup.hpp"
// #include "matplotlibcpp.h"

#define UNIT_SECONDS microseconds
#define LABEL_SECONDS "μs"

#define START_TIMER() begin = chrono::steady_clock::now()
#define END_TIMER() end = chrono::steady_clock::now()
#define ADD_TIME_TO(x) x += chrono::duration_cast<chrono::UNIT_SECONDS>(end - begin).count()

#define EPS 1e-5

using namespace std;

double sample(double x, double y) {
  if (x < 0) x += 2 * M_PI;  // periodic extension
  if (y < 0) y += 2 * M_PI;  // periodic extension
  return sin(98 * x);
  //  + cos(49 * y + 32 * x) + cos(71 * x + 17 * y);
}

double diff_x_sample(double x, double y) {
  if (x < 0) x += 2 * M_PI;  // periodic extension
  if (y < 0) y += 2 * M_PI;  // periodic extension
  return 98 * cos(98 * x);
  //  - 32 * sin(49 * y + 32 * x) - 71 * sin(71 * x + 17 * y);
}

double diff_y_sample(double x, double y) {
  if (x < 0) x += 2 * M_PI;  // periodic extension
  if (y < 0) y += 2 * M_PI;  // periodic extension
  return -49 * sin(49 * y + 32 * x) - 17 * sin(71 * x + 17 * y);
}

int main(void) {
  const string filename = "data/data.txt";   // name of the output file
  const int nx = 16, ny = 16;                // number of points in x and y (must be a power of 2)
  const int nu1 = 1, nu2 = 1;                // parameters of the system
  const int dim = 2 * nx * ny;               // dimension of the system (the '2' is for the real and imaginary part)
  double t = 0.;                             // initial time of integration
  const double T = 5.;                       // final time of integration
  const double h = 0.005;                    // initial step size
  const int maxNumSteps = (int)(T / h) + 1;  // maximum number of steps
  int numSteps = 0;                          // number of steps
  Vector<double> x(dim);                     // initial condition
  Vector<double> aux(dim);

  chrono::steady_clock::time_point begin, end;
  int64_t total_write = 0, total_computations = 0;
  // begin_computations = chrono::steady_clock::now();

  START_TIMER();
  Args prm(nx, ny, h, nu1, nu2);  // parameters for the system (nu1, nu2, nn, k1, k2, tmp
  set_data(x, prm);
  END_TIMER();
  ADD_TIME_TO(total_computations);

  // file handling
  START_TIMER();
  ofstream file;  // output file
  file.open(filename);
  write(x, prm.nn, t, file);
  END_TIMER();
  ADD_TIME_TO(total_write);

  Vector<double> diff_y(dim), diff_y_aux(dim), tmp(dim);
  double a1 = 2 * M_PI / prm.nx;
  double a2 = 2 * M_PI / prm.ny;
  for (int i = 0; i < prm.nx; i++) {
    for (int j = 0; j < 2 * prm.ny; j += 2) {
      x[2 * i * prm.ny + j] = sample(prm.k1[2 * i * prm.ny + j] * a1, prm.k2[2 * i * prm.ny + j] * a2);
      diff_y[2 * i * prm.ny + j] = diff_x_sample(prm.k1[2 * i * prm.ny + j] * a1, prm.k2[2 * i * prm.ny + j] * a2);
      x[2 * i * prm.ny + j + 1] = 0;
      diff_y[2 * i * prm.ny + j + 1] = 0;
    }
  }

  // compute the FFT of x
  START_TIMER();
  fft_n(x, prm.nn, 1, 1);
  END_TIMER();
  ADD_TIME_TO(total_computations);

  // // we want to compare different methods to compute the transform of (u_x)^2
  // int repetitions = 100, count = 0;
  // int64_t total_real = 0, total_2 = 0, total_3 = 0;

  // // while (count < repetitions) {
  // // 1. Real solution F((u_x)^2)

  // file << "k1" << endl;
  // write(prm.k1, prm.nn, t, file, true);

  // file << "F(u_y)" << endl;
  // write(diff_y, prm.nn, t, file, true);
  // fft_n(diff_y, prm.nn, 1, 1);
  // write(diff_y, prm.nn, t, file, true);

  // diff_y_aux = diff_y;
  // START_TIMER();
  // diff_y_aux *= diff_y_aux;
  // fft_n(diff_y_aux, prm.nn, 1, 1);
  // END_TIMER();
  // ADD_TIME_TO(total_real);

  // // 2. F([F⁻¹(i * k₁ * û)]²)
  // START_TIMER();
  // tmp = prm.k1 * x;
  // double tmp1;
  // for (int i = 0; i < dim; i += 2) {  // multiply by i
  //   tmp1 = tmp[i];
  //   tmp[i] = -tmp[i + 1];
  //   tmp[i + 1] = tmp1;
  // }
  // file << "F(u_y) approx" << endl;
  // write(tmp, prm.nn, t, file, true);
  // fft_n(tmp, prm.nn, -1, 1);
  // write(tmp, prm.nn, t, file, true);
  // tmp *= tmp;
  // fft_n(tmp, prm.nn, 1, 1);
  // END_TIMER();
  // ADD_TIME_TO(total_2);

  // 3. F((u_x)^2) = F(u_x) * F(u_x) = - (k1 . û) * (k1 . û)
  // START_TIMER();
  // tmp = prm.k1 * x;
  // // do the convolution

  // for (int i = 0; i < nx; i++){
  //   for (int j = 0; j < 2 * ny; j
  // }
  //   count++;
  // }

  // file << "diff (fourier_space): " << endl;
  // diff_y -= tmp;
  // write(diff_y_aux, prm.nn, t, file, true);

  do {
    START_TIMER();
    stepFiniteDiff(x, t, h, prm);
    aux = x;                    // copy x to aux
    fft_n(aux, prm.nn, -1, 1);  // inverse FFT to output the solution
    END_TIMER();
    ADD_TIME_TO(total_computations);

    START_TIMER();
    write(aux, prm.nn, t, file);  // write the solution into the file
    END_TIMER();
    ADD_TIME_TO(total_write);

    numSteps++;
  } while (t < T - EPS && numSteps <= maxNumSteps);

  cout << "Total time for computations: " << total_computations << " " << LABEL_SECONDS << endl;
  cout << "Total time for writing: " << total_write << " " << LABEL_SECONDS << endl;

  file.close();
  return 0;
}