#include <fftw3.h>

#include <chrono>
#include <cstring>
#include <fstream>
#include <iostream>

#include "../include/energy.hpp"
#include "../include/misc.hpp"
#include "../include/scheme_fd.hpp"
#include "../include/vector.hpp"

#define UNIT_SECONDS 1000  // in microseconds
#define LABEL_SECONDS "ms"

#define START_TIMER() begin = chrono::steady_clock::now()
#define END_TIMER() end = chrono::steady_clock::now()
#define ADD_TIME_TO(x) x += chrono::duration_cast<chrono::microseconds>(end - begin).count()
#define FFTW_FLAG FFTW_ESTIMATE  // options: FFTW_ESTIMATE, FFTW_MEASURE.
// FFTW_MEASURE instructs FFTW to run and measure the execution time of several FFTs in order to find the best way to compute the transform of size n.This process takes some time(usually a few seconds), depending on your machine and on the size of the transform.FFTW_ESTIMATE, on the contrary, does not run any computation and just builds a reasonable plan that is probably sub - optimal.In short, if your program performs many transforms of the same size and initialization time is not important, use FFTW_MEASURE; otherwise use the estimate.
#define PLOT() ((write_solution) && (numSteps % freq_plot_sol == 0))
#define WRITE() ((write_energy) && (numSteps % freq_plot_e == 0))

#define EPS 1e-5
using namespace std;

int main(void) {
  // ----------------- PARAMETERS ------------------
  const string filename_solution = "data/solution.txt";  // name of the output file
  const string filename_energy = "data/energy.txt";      // name of the output file
  const string filename_input = "data/input.txt";        // name of the input file
  const string space = "    ";                           // space to print
  uint count = 5;                                        // progress percentage interval to print (each count%)
  chrono::steady_clock::time_point begin, end;           // variables to measure the time
  int64_t total_write = 0, total_computations = 0, total_fft_precomputation = 0;
  int nx, ny;       // number of points in x and y (must be a power of 2), do NOT use uint fftw uses int
  double T;         // final time of integration
  double dt;        // initial step size
  double nu1, nu2;  // parameters of the system (by default)
  // DO NOT DECREASE THE NU'S BELOW 0.2 BECAUSE THE SOLUTION WILL START TO EXPLODE
  bool averaged_solution;  // 1 if we want to average the solution, 0 otherwise
  bool write_energy;       // whether to write the energy to a file or not
  bool write_solution;     // whether to write the solution to a file or not
  uint freq_plot_sol;      // frequency of the plot of the solution
  uint freq_plot_e;        // frequency of the plot of the energy
  // -----------------------------------------------

  START_TIMER();
  // ------------- File input setup ----------------
  ifstream file_input;
  file_input.open("data/input.txt");
  string tmp;
  if (file_input.is_open()) {
    file_input >> tmp >> nx;
    file_input >> tmp >> ny;
    file_input >> tmp >> nu1;
    file_input >> tmp >> nu2;
    file_input >> tmp >> T;
    file_input >> tmp >> dt;
    file_input >> tmp >> write_solution;
    file_input >> tmp >> write_energy;
    file_input >> tmp >> averaged_solution;
    file_input >> tmp >> freq_plot_sol;
    file_input >> tmp >> freq_plot_e;
  }
  file_input.close();
  // -----------------------------------------------

  // plot the parameters
  cout << "Size of the grid:      " << space << nx << " x " << ny << endl;
  cout << "Parameters:            " << space << "nu1 = " << nu1 << endl;
  cout << "                       " << space << "nu2 = " << nu2 << endl;
  cout << "Final time:            " << space << T << endl;
  cout << "Step size:             " << space << dt << endl;
  cout << "Use averaged solution? " << space << (averaged_solution ? "yes" : "no") << endl;
  cout << "Write solution to file?" << space << (write_solution ? "yes" : "no") << endl;
  cout << "Write energy to file?  " << space << (write_energy ? "yes" : "no") << endl;
  END_TIMER();
  ADD_TIME_TO(total_write);

  // ----------------- PARAMETERS ------------------
  int ny_complex = ny / 2 + 1;                  // number of points in y in the frequency space to optimize the memory in the FFT
  uint dim = (uint)(nx * ny);                   // dimension of the physical system
  uint dim_f = (uint)(nx * ny_complex);         // dimension of the system in the frequency space
  double fraction_completed = T / 100.;         // fraction of the integration time to print
  const uint maxNumSteps = (uint)(T / dt) + 1;  // maximum number of steps (I add 1 to be sure to have at least T/dt steps)
  double t = 0.;                                // initial time of integration
  uint numSteps = 0;                            // number of steps
  // vector to store all the solutions
  double energy, energy2;
  double *kx = (double *)malloc(dim_f * sizeof(double));
  double *ky = (double *)malloc(dim_f * sizeof(double));
  double *C = (double *)malloc(dim_f * sizeof(double));
  // -----------------------------------------------

  // ----------------- FFT SETUP -------------------
  double *u = fftw_alloc_real(dim);  // u is the solution,
  double *aux_x = fftw_alloc_real(dim);
  double *aux_y = fftw_alloc_real(dim);
  fftw_complex *Fu = fftw_alloc_complex(dim_f);  // Fu is the DFT of u
  fftw_complex *Fu_aux = fftw_alloc_complex(dim_f);
  fftw_complex *aux_x_fourier = fftw_alloc_complex(dim_f);
  fftw_complex *aux_y_fourier = fftw_alloc_complex(dim_f);

  fftw_plan p_for_u, p_back_u, p_for_aux_x, p_for_aux_y, p_back_aux_x, p_back_aux_y;

  START_TIMER();
  // Remember that for inverse transforms, that is for _c2r_ transforms, the input array is overwritten, so you should copy it before if you want to preserve the input data.
  // Also remember that the fftw computes unnormalized transforms, so doing a forward transform followed by a backward transform will multiply the input by n.
  p_for_u = fftw_plan_dft_r2c_2d(nx, ny, u, Fu, FFTW_FLAG);                      // forward
  p_back_u = fftw_plan_dft_c2r_2d(nx, ny, Fu_aux, u, FFTW_FLAG);                 // backward
  p_for_aux_x = fftw_plan_dft_r2c_2d(nx, ny, aux_x, aux_x_fourier, FFTW_FLAG);   // forward
  p_for_aux_y = fftw_plan_dft_r2c_2d(nx, ny, aux_y, aux_y_fourier, FFTW_FLAG);   // forward
  p_back_aux_x = fftw_plan_dft_c2r_2d(nx, ny, aux_x_fourier, aux_x, FFTW_FLAG);  // backward
  p_back_aux_y = fftw_plan_dft_c2r_2d(nx, ny, aux_y_fourier, aux_y, FFTW_FLAG);  // backward
  END_TIMER();
  ADD_TIME_TO(total_fft_precomputation);
  // -----------------------------------------------

  START_TIMER();
  set_data(u, nx, ny);
  set_wave_numbers(kx, ky, nx, ny_complex);
  set_C(C, kx, ky, nx, ny_complex, dt, nu1, nu2);
  fftw_execute(p_for_u);  // Now Fu contains the DFT of u
  energy = E(u, nx, ny);
  END_TIMER();
  ADD_TIME_TO(total_computations);

  // file handling
  START_TIMER();
  ofstream file_sol, file_E;                       // output files
  if (averaged_solution && (PLOT() || WRITE())) {  // Fu[0], the fourier coefficient with k_1 = k_2 = 0, is the mean of the solution
    for (uint i = 0; i < dim; i++) u[i] -= Fu[0][0];
  }
  if (PLOT()) {
    file_sol.open(filename_solution);
    write(u, nx, ny, t, file_sol);  // write the solution into the file
  }
  if (WRITE()) {  // write the energy into the file
    file_E.open(filename_energy);
    // in the 4th column plot u(pi, pi)
    file_E << t << " " << energy << " 0.0 " << u[(ny / 2) * (nx + 1)] << endl;  // write the energy into the file
  }
  END_TIMER();
  ADD_TIME_TO(total_write);
  numSteps++;
  do {
    START_TIMER();
    // stepIMEXRK4(x, t, prm);
    // stepFiniteDiff(x, t, prm);
    {  // step in finite differences
       // multiply Fu by i * k1 and store the result in aux_x_fourier and multiply Fu by i * k2 and store the result in aux_y_fourier
      for (uint i = 0; i < dim_f; i++) {
        aux_x_fourier[i][0] = -Fu[i][1] * kx[i];
        aux_x_fourier[i][1] = Fu[i][0] * kx[i];
        aux_y_fourier[i][0] = -Fu[i][1] * ky[i];
        aux_y_fourier[i][1] = Fu[i][0] * ky[i];
      }
      // here we lose the information of aux_x_fourier and aux_y_fourier, but we don't need it anymore (we will recover it later)
      fftw_execute(p_back_aux_x);  // Now aux_x contains the inverse DFT of aux_x_fourier, which is u_x
      fftw_execute(p_back_aux_y);  // Now aux_y contains the inverse DFT of aux_y_fourier, which is u_y
      // normalize the inverse DFT
      for (uint i = 0; i < dim; i++) {
        aux_x[i] /= dim;
        aux_y[i] /= dim;
      }

      // compute (u_x)^2 and (u_y)^2
      // if (numSteps == 0)
      //   write(aux_x, nx, ny, t, file_E);
      for (uint i = 0; i < dim; i++) {
        aux_x[i] *= aux_x[i];
        aux_y[i] *= aux_y[i];
      }
      // compute the DFT of (u_x)^2 and store the result in aux_x_fourier and compute the DFT of (u_y)^2 and store the result in aux_y_fourier (here we recover the information of aux_x_fourier and aux_y_fourier)
      fftw_execute(p_for_aux_x);
      fftw_execute(p_for_aux_y);
      // if (numSteps == 0)
      //   write_old(Fu, nx, ny_complex, t, file_E);
      // compute the new iterate of the finite differences scheme
      for (uint i = 0; i < dim_f; i++) {
        Fu[i][0] = C[i] * (Fu[i][0] - dt * 0.5 * (aux_x_fourier[i][0] + nu2 / nu1 * aux_y_fourier[i][0]));
        Fu[i][1] = C[i] * (Fu[i][1] - dt * 0.5 * (aux_x_fourier[i][1] + nu2 / nu1 * aux_y_fourier[i][1]));
      }
      // increase the time
      t += dt;
    }
    // normalize the inverse DFT. We first substract the mean of the solution and then we normalize it
    // for (uint i = 0; i < dim; i++) u[i] /= dim;
    if (PLOT() || WRITE()) {
      // we would lose the information of Fu, so we save it in Fu_aux and now we will lose the information of Fu_aux
      memcpy(Fu_aux, Fu, dim_f * sizeof(fftw_complex));
      fftw_execute(p_back_u);                                           // Now u contains the inverse DFT of Fu, which is the new iterate of the solution in the physical space
      if (averaged_solution) {                                          // Fu[0], the fourier coefficient with k_1 = k_2 = 0, is the mean of the solution
        for (uint i = 0; i < dim; i++) u[i] = (u[i] - Fu[0][0]) / dim;  // we need to normalize Fu accordingly in order to be the mean of the solution
      } else {
        for (uint i = 0; i < dim; i++) u[i] /= dim;  // we need to normalize Fu accordingly in order to be the mean of the solution
      }
    }
    if (WRITE()) {
      energy2 = energy;
      energy = E(u, nx, ny);
    }
    END_TIMER();
    ADD_TIME_TO(total_computations);

    START_TIMER();
    if (PLOT())
      write(u, nx, ny, t, file_sol);  // write the solution into the file
    if (WRITE())                      // write the energy into the file
      file_E << t << " " << energy << " " << (energy - energy2) / dt << " " << u[(ny / 2) * (nx + 1)] << endl;
    if (t > fraction_completed * count - EPS) {
      cout << count << "%" << endl;
      count += 5;
    }
    END_TIMER();
    ADD_TIME_TO(total_write);

    numSteps++;
  } while (t < T - EPS && numSteps <= maxNumSteps + 1);  // We add 1 because we count the initial condition as a step

  // creation of tmp files to comunicate with bash scripts in order to plot or not the solution and the energy
  START_TIMER();
  ofstream file_tmp;
  if (write_solution) {
    file_tmp.open("data/tmp_write_sol.txt");
    file_tmp.close();
  }
  if (write_energy) {
    file_tmp.open("data/tmp_write_E.txt");
    file_tmp.close();
  }

  cout << "Total time for fftw precomputation: " << total_fft_precomputation / UNIT_SECONDS << " " << LABEL_SECONDS << endl;
  cout << "Total time for computations: " << total_computations / UNIT_SECONDS << " " << LABEL_SECONDS << endl;
  cout << "Total time for read/write: " << total_write / UNIT_SECONDS << " " << LABEL_SECONDS << endl;

  if (write_solution) file_sol.close();
  if (write_energy) file_E.close();

  // clean variables
  free(kx);
  free(ky);
  free(C);

  // -----------------------------------------------
  // ----------------- FFT CLEANUP -----------------
  // -----------------------------------------------
  fftw_destroy_plan(p_for_u);
  fftw_destroy_plan(p_back_u);
  fftw_destroy_plan(p_for_aux_x);
  fftw_destroy_plan(p_for_aux_y);
  fftw_destroy_plan(p_back_aux_x);
  fftw_destroy_plan(p_back_aux_y);
  fftw_free(u);
  fftw_free(aux_x);
  fftw_free(aux_y);
  fftw_free(Fu);
  fftw_free(Fu_aux);
  fftw_free(aux_x_fourier);
  fftw_free(aux_y_fourier);
  fftw_cleanup();
  // -----------------------------------------------
  // -----------------------------------------------
  // -----------------------------------------------
  return 0;
}