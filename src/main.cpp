#include <fftw3.h>
#include <python3.11/Python.h>

#include <chrono>
#include <fstream>
#include <iostream>

#include "../include/energy.hpp"
#include "../include/fft.hpp"
#include "../include/misc.hpp"
#include "../include/scheme_fd.hpp"
#include "../include/vector.hpp"

#define UNIT_SECONDS 1000  // in microseconds
#define LABEL_SECONDS "ms"

#define START_TIMER() begin = chrono::steady_clock::now()
#define END_TIMER() end = chrono::steady_clock::now()
#define ADD_TIME_TO(x) x += chrono::duration_cast<chrono::microseconds>(end - begin).count()

#define EPS 1e-5
using namespace std;

int main(int argc, char const *argv[]) {
  // -----------------------------------------------
  // ----------------- PARAMETERS ------------------
  // -----------------------------------------------
  const bool averaged_solution = true;               // 1 if we want to average the solution, 0 otherwise
  const string filename_solution = "data/data.txt";  // name of the output file
  const string filename_energy = "data/energy.txt";  // name of the output file
  const uint nx = 32, ny = 32;                       // number of points in x and y (must be a power of 2)
  const double T = 500.;                             // final time of integration
  const double fraction_completed = T / 100.;        // fraction of the integration time to print
  uint count = 5;
  const double dt = 0.001;            // initial step size
  double nu1 = 0.5, nu2 = 0.5;        // parameters of the system (by default)
  const bool write_solution = false;  // whether to write the solution to a file or not
  const bool write_energy = true;     // whether to write the energy to a file or not
  // DO NOT DECREASE THE NU'S BELOW 0.2 BECAUSE THE SOLUTION WILL START TO EXPLODE
  if (argc > 1) nu1 = atof(argv[1]);  // if the user specifies the parameters, we use them
  if (argc > 2) nu2 = atof(argv[2]);  // if the user specifies the parameters, we use them
  // -----------------------------------------------
  // -----------------------------------------------
  // -----------------------------------------------

  // -----------------------------------------------
  // ----------------- FFT SETUP -------------------
  // -----------------------------------------------
  // fftw_complex *in, *out;
  // fftw_plan p;

  // in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * nx * ny);
  // out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * nx * ny);

  // p = fftw_plan_dft_2d(nx, ny, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

  const uint dim = 2 * nx * ny;                // dimension of the system (the '2' is for the real and imaginary part)
  double t = 0.;                               // initial time of integration
  const uint maxNumSteps = (int)(T / dt) + 1;  // maximum number of steps (I add 1 to be sure to have at least T/dt steps)
  uint numSteps = 0;                           // number of steps
  vector<double> x(dim);                       // initial condition
  vector<double> aux(dim);
  // vector to store all the solutions
  double energy, energy2;

  chrono::steady_clock::time_point begin, end;
  int64_t total_write = 0, total_computations = 0;
  // begin_computations = chrono::steady_clock::now();

  START_TIMER();
  Args prm(nx, ny, dt, nu1, nu2);  // parameters for the system (nu1, nu2, nn, k1, k2, tmp
  set_data(x, prm);
  aux = x;
  fft_n(x, prm.nn, 1, 1);              // compute the FFT of x
  if (averaged_solution) aux -= x[0];  // x[0], the fourier coefficient with k_1 = k_2 = 0, is the mean of the solution
  energy = E(aux, prm);
  END_TIMER();
  ADD_TIME_TO(total_computations);

  // file handling
  START_TIMER();
  ofstream file_sol, file_E;  // output files
  file_sol.open(filename_solution);
  file_E.open(filename_energy);
  if (write_solution)
    write(aux, prm.nn, t, file_sol);
  if (write_energy)
    file_E << t << " " << energy << " 0.0" << endl;  // write the energy into the file
  END_TIMER();
  ADD_TIME_TO(total_write);

  do {
    START_TIMER();
    // stepIMEXRK4(x, t, prm);
    stepFiniteDiff(x, t, prm);
    aux = x;                             // copy x to aux
    fft_n(aux, prm.nn, -1, 1);           // inverse FFT to output the solution
    if (averaged_solution) aux -= x[0];  // x[0], the fourier coefficient with k_1 = k_2 = 0, is the mean of the solution
    energy2 = energy;
    energy = E(aux, prm);
    END_TIMER();
    ADD_TIME_TO(total_computations);

    START_TIMER();
    if (write_solution)
      write(aux, prm.nn, t, file_sol);  // write the solution into the file
    if (write_energy)
      file_E << t << " " << energy << " " << (energy - energy2) / dt << endl;  // write the energy into the file
    if (t > fraction_completed * count) {
      cout << "Completed " << count << "% of the integration" << endl;
      count += 5;
    }
    END_TIMER();
    ADD_TIME_TO(total_write);

    numSteps++;
  } while (t < T - EPS && numSteps <= maxNumSteps);

  cout << "Total time for computations: " << total_computations / UNIT_SECONDS << " " << LABEL_SECONDS << endl;
  cout << "Total time for writing: " << total_write / UNIT_SECONDS << " " << LABEL_SECONDS << endl;

  file_sol.close();
  file_E.close();
  return 0;
}