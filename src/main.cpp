#include <chrono>
#include <fstream>
#include <iostream>

#include "../include/fft.hpp"
#include "../include/misc.hpp"
#include "../include/vector.hpp"

#define UNIT_SECONDS microseconds
#define LABEL_SECONDS "μs"

#define START_TIMER() begin = chrono::steady_clock::now()
#define END_TIMER() end = chrono::steady_clock::now()
#define ADD_TIME_TO(x) x += chrono::duration_cast<chrono::UNIT_SECONDS>(end - begin).count()

#define EPS 1e-5
using namespace std;

int main(void) {
  // -----------------------------------------------
  // ----------------- PARAMETERS ------------------
  // -----------------------------------------------
  const bool averaged_solution = true;      // 1 if we want to average the solution, 0 otherwise
  const string filename = "data/data.txt";  // name of the output file
  const uint nx = 32, ny = 32;              // number of points in x and y (must be a power of 2)
  const double nu1 = 0.05, nu2 = 0.05;      // parameters of the system
  const double T = 5.;                      // final time of integration
  const double h = 0.005;                   // initial step size
  // -----------------------------------------------
  // -----------------------------------------------
  // -----------------------------------------------

  const uint dim = 2 * nx * ny;               // dimension of the system (the '2' is for the real and imaginary part)
  double t = 0.;                              // initial time of integration
  const uint maxNumSteps = (int)(T / h) + 1;  // maximum number of steps
  uint numSteps = 0;                          // number of steps
  vector<double> x(dim);                      // initial condition
  vector<double> aux(dim);

  chrono::steady_clock::time_point begin, end;
  int64_t total_write = 0, total_computations = 0;
  // begin_computations = chrono::steady_clock::now();

  START_TIMER();
  Args prm(nx, ny, h, nu1, nu2);  // parameters for the system (nu1, nu2, nn, k1, k2, tmp
  set_data(x, prm);
  if (averaged_solution) x -= mean(x, prm);
  END_TIMER();
  ADD_TIME_TO(total_computations);

  // file handling
  START_TIMER();
  ofstream file;  // output file
  file.open(filename);
  write(x, prm.nn, t, file);
  END_TIMER();
  ADD_TIME_TO(total_write);

  // compute the FFT of x
  START_TIMER();
  fft_n(x, prm.nn, 1, 1);
  END_TIMER();
  ADD_TIME_TO(total_computations);

  do {
    START_TIMER();
    stepFiniteDiff(x, t, prm);
    aux = x;                    // copy x to aux
    fft_n(aux, prm.nn, -1, 1);  // inverse FFT to output the solution
    if (averaged_solution) aux -= mean(aux, prm);
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