#include <fftw3.h>

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
#define FFTW_FLAG FFTW_ESTIMATE  // other options: FFTW_MEASURE. TFFTW_MEASURE
// FFTW_MEASURE instructs FFTW to run and measure the execution time of several FFTs in order to find the best way to compute the transform of size n.This process takes some time(usually a few seconds), depending on your machine and on the size of the transform.FFTW_ESTIMATE, on the contrary, does not run any computation and just builds a reasonable plan that is probably sub - optimal.In short, if your program performs many transforms of the same size and initialization time is not important, use FFTW_MEASURE; otherwise use the estimate.
#define EPS 1e-5
using namespace std;

// #define COPY(x, y)
//   for (uint i = 0; i < nx * ny; i++) {
//     (y)[i][0] = (x)[i][0];
//     (y)[i][1] = (x)[i][1];
//   }

int main(int argc, char const *argv[]) {
  // -----------------------------------------------
  // ----------------- PARAMETERS ------------------
  // -----------------------------------------------
  const bool averaged_solution = true;               // 1 if we want to average the solution, 0 otherwise
  const string filename_solution = "data/data.txt";  // name of the output file
  const string filename_energy = "data/energy.txt";  // name of the output file
  const uint nx = 16, ny = 16;                       // number of points in x and y (must be a power of 2)
  const uint ny_complex = ny / 2 + 1;                // number of points in y in the frequency space to optimize the memory in the FFT
  const double T = 500.;                             // final time of integration
  const double fraction_completed = T / 100.;        // fraction of the integration time to print
  uint count = 5;
  const double dt = 0.001;           // initial step size
  double nu1 = 0.5, nu2 = 0.5;       // parameters of the system (by default)
  const bool write_solution = true;  // whether to write the solution to a file or not
  const bool write_energy = true;    // whether to write the energy to a file or not
  // DO NOT DECREASE THE NU'S BELOW 0.2 BECAUSE THE SOLUTION WILL START TO EXPLODE
  if (argc > 1) nu1 = atof(argv[1]);  // if the user specifies the parameters, we use them
  if (argc > 2) nu2 = atof(argv[2]);  // if the user specifies the parameters, we use them
  // -----------------------------------------------
  // -----------------------------------------------
  // -----------------------------------------------
  double t = 0.;                               // initial time of integration
  const uint maxNumSteps = (int)(T / dt) + 1;  // maximum number of steps (I add 1 to be sure to have at least T/dt steps)
  uint numSteps = 0;                           // number of steps
  // vector to store all the solutions
  double energy, energy2;
  chrono::steady_clock::time_point begin, end;
  int64_t total_write = 0, total_computations = 0;

  // -----------------------------------------------
  // ----------------- FFT SETUP -------------------
  // -----------------------------------------------
  double *u = fftw_alloc_real(nx * ny);  // u is the solution,
  double *aux_x = fftw_alloc_real(nx * ny);
  double *aux_y = fftw_alloc_real(nx * ny);
  fftw_complex *Fu = fftw_alloc_complex(nx * ny_complex);  // Fu is the DFT of u
  fftw_complex *aux_x_fourier = fftw_alloc_complex(nx * ny_complex);
  fftw_complex *aux_y_fourier = fftw_alloc_complex(nx * ny_complex);

  fftw_plan p_for_u, p_back_u, p_for_aux_x, p_for_aux_y, p_back_aux_x, p_back_aux_y;

  START_TIMER();
  // Remember that for inverse transforms, that is for _c2r_ transforms, the input array is overwritten, so you should copy it before if you want to preserve the input data.
  p_for_u = fftw_plan_dft_r2c_2d(nx, ny, u, Fu, FFTW_FLAG);
  p_back_u = fftw_plan_dft_c2r_2d(nx, ny, Fu, u, FFTW_FLAG);
  p_for_aux_x = fftw_plan_dft_r2c_2d(nx, ny, aux_x, aux_x_fourier, FFTW_FLAG);
  p_for_aux_y = fftw_plan_dft_r2c_2d(nx, ny, aux_y, aux_y_fourier, FFTW_FLAG);
  p_back_aux_x = fftw_plan_dft_c2r_2d(nx, ny, aux_x_fourier, aux_x, FFTW_FLAG);
  p_back_aux_y = fftw_plan_dft_c2r_2d(nx, ny, aux_y_fourier, aux_y, FFTW_FLAG);
  // -----------------------------------------------
  // -----------------------------------------------
  // -----------------------------------------------

  // Args prm(nx, ny, dt, nu1, nu2);  // parameters for the system (nu1, nu2, nn, k1, k2, tmp
  set_data(u, nx, ny);
  double *k1 = (double *)malloc(nx * ny_complex * sizeof(double));
  double *k2 = (double *)malloc(nx * ny_complex * sizeof(double));
  set_wave_numbers(k1, k2, nx, ny_complex);
  double a1 = 2 * M_PI / nx;
  double a2 = 2 * M_PI / ny;
  for (uint i = 0; i < nx; i++) {
    for (uint j = 0; j < ny; j++) {
      aux_y[i * ny + j] = d_y_u0(i * a1, j * a2);
    }
  }
  // if (averaged_solution) {  // Fu[0], the fourier coefficient with k_1 = k_2 = 0, is the mean of the solution
  //   for (uint i = 0; i < nx * ny; i++) u[i] -= Fu[0][0];
  // }
  energy = E(u, nx, ny);
  END_TIMER();
  ADD_TIME_TO(total_computations);
  ofstream file_sol, file_E;  // output files
  file_sol.open(filename_solution);
  file_E.open(filename_energy);
  file_sol << "nx = " << nx << endl;
  file_sol << "ny = " << ny << endl;
  file_sol << "ny_complex = " << ny_complex << endl;
  file_sol << "u = (before fftw_execute call)" << endl;
  write(u, nx, ny, t, file_sol);
  // file_sol << "u before fftw_execute call" << endl;
  // write(u, nx, ny, t, file_sol);
  fftw_execute(p_for_u);  // Now Fu contains the DFT of u
  // file_sol << "u after fftw_execute call" << endl;
  // write(u, nx, ny, t, file_sol);
  // file_sol << "Fu before fftw_execute call" << endl;
  // write_old(Fu, nx, ny_complex, t, file_sol);
  // fftw_execute(p_back_u);  // Now u contains the inverse DFT of Fu, which is the new iterate of the solution in the physical space
  // file_sol << "Fu after fftw_execute call" << endl;
  // write_old(Fu, nx, ny_complex, t, file_sol);
  // file_sol << "u after fftw_execute call" << endl;
  // write(u, nx, ny, t, file_sol);
  // normalization of u
  // for (uint i = 0; i < nx * ny; i++) {
  //   u[i] /= nx * ny;
  // }
  // file_sol << "u after normalization" << endl;
  // write(u, nx, ny, t, file_sol);

  // file_sol << "u after fftw_execute call" << endl;
  // write(u, nx, ny, t, file_sol);
  for (uint i = 0; i < nx * ny_complex; i++) {
    aux_x_fourier[i][0] = -Fu[i][1] * k2[i];
    aux_x_fourier[i][1] = Fu[i][0] * k2[i];
  }
  vector<double> U(2 * nx * ny);
  for (uint i = 0; i < nx * ny; i++) {
    U[2 * i] = u[i];
    U[2 * i + 1] = 0;
  }
  vector<uint> nn = {nx, ny};
  fft_n(U, nn, 1, -1);
  START_TIMER();
  file_sol << "fft(u) with old style" << endl;
  write_old_old(U, nn, t, file_sol, true);
#define round_to_zero(x) (fabs(x) < 1e-10 ? 0 : x)
  file_sol << "fft(u) with new style" << endl;
  write_old(Fu, nx, ny_complex, t, file_sol);
  file_sol << "fft(u_y) with new style (using k2) (before)" << endl;
  write_old(aux_x_fourier, nx, ny_complex, t, file_sol);
  for (uint i = 0; i < nx * ny; i++) {
    U[2 * i] = aux_y[i];
    U[2 * i + 1] = 0;
  }
  fft_n(U, nn, 1, -1);
  file_sol << "fft(u_y) with new style (using u_y)" << endl;
  fftw_execute(p_for_aux_y);
  write_old(aux_y_fourier, nx, ny_complex, t, file_sol);
  file_sol << "fft(u_y) with old style" << endl;
  write_old_old(U, nn, t, file_sol, true);
  fftw_execute(p_back_aux_x);  // Now aux_x contains the inverse DFT of aux_x_fourier, which is u_x
  // file_sol << "fft(u_x) with new style (after)" << endl;
  // write_old(aux_x_fourier, nx, ny_complex, t, file_sol);
  // file handling
  file_sol << "u_x = " << endl;
  write(aux_y, nx, ny, t, file_sol);
  file_sol << "u_x with new style after fftw_execute call" << endl;
  write(aux_x, nx, ny, t, file_sol);
  file_sol << "difference of u_x (should be = 0)" << endl;
  double *diff = (double *)malloc(nx * ny * sizeof(double));
  for (uint i = 0; i < nx * ny; i++) {
    diff[i] = round_to_zero(aux_x[i] / (nx * ny) - aux_y[i]);
  }
  write(diff, nx, ny, t, file_sol);
  // if (write_energy)
  //   file_E << t << " " << energy << " 0.0" << endl;  // write the energy into the file
  END_TIMER();
  ADD_TIME_TO(total_write);

  // cout << "k_x = " << k1 << endl;
  // cout << "k_y = " << k2 << endl;
#undef round_to_zero

  // do {
  //   START_TIMER();
  //   // stepIMEXRK4(x, t, prm);
  //   // stepFiniteDiff(x, t, prm);
  //   {  // step in finite differences
  //      // multiply Fu by i * k1 and store the result in aux_x_fourier and multiply Fu by i * k2 and store the result in aux_y_fourier
  //     for (uint i = 0; i < nx * ny_complex; i++) {
  //       aux_x_fourier[i][0] = -Fu[i][1] * prm.k1[i];
  //       aux_x_fourier[i][1] = Fu[i][0] * prm.k1[i];
  //       aux_y_fourier[i][0] = -Fu[i][1] * prm.k2[i];
  //       aux_y_fourier[i][1] = Fu[i][0] * prm.k2[i];
  //     }
  //     fftw_execute(p_back_aux_x);  // Now aux_x contains the inverse DFT of aux_x_fourier, which is u_x
  //     fftw_execute(p_back_aux_y);  // Now aux_y contains the inverse DFT of aux_y_fourier, which is u_y
  //     // compute (u_x)^2 and (u_y)^2
  //     for (uint i = 0; i < nx * ny; i++) {
  //       aux_x[i] *= aux_x[i];
  //       aux_y[i] *= aux_y[i];
  //     }
  //     // compute the DFT of (u_x)^2 and store the result in aux_x_fourier and compute the DFT of (u_y)^2 and store the result in aux_y_fourier
  //     fftw_execute(p_for_aux_x);
  //     fftw_execute(p_for_aux_y);
  //     // compute the new iterate of the finite differences scheme
  //     for (uint i = 0; i < nx * ny_complex; i++) {
  //       Fu[i][0] = prm.tmp[i] * (Fu[i][0] - prm.h * (0.5 * (aux_x_fourier[i][0] + ((prm.nu2 / prm.nu1) * aux_y_fourier[i][0]))));
  //       Fu[i][1] = prm.tmp[i] * (Fu[i][1] - prm.h * (0.5 * (aux_x_fourier[i][1] + ((prm.nu2 / prm.nu1) * aux_y_fourier[i][1]))));
  //     }
  //     // increase the time
  //     t += prm.h;
  //   }
  //   fftw_execute(p_back_u);   // Now u contains the inverse DFT of Fu, which is the new iterate of the solution in the physical space
  //   if (averaged_solution) {  // Fu[0], the fourier coefficient with k_1 = k_2 = 0, is the mean of the solution
  //     for (uint i = 0; i < nx * ny; i++) u[i] -= Fu[0][0];
  //   }
  //   energy2 = energy;
  //   energy = E(u, prm);
  //   END_TIMER();
  //   ADD_TIME_TO(total_computations);

  //   START_TIMER();
  //   if (write_solution)
  //     write(u, prm.nn, t, file_sol);  // write the solution into the file
  //   if (write_energy)
  //     file_E << t << " " << energy << " " << (energy - energy2) / dt << endl;  // write the energy into the file
  //   if (t > fraction_completed * count) {
  //     cout << "Completed " << count << "% of the integration" << endl;
  //     count += 5;
  //   }
  //   END_TIMER();
  //   ADD_TIME_TO(total_write);

  //   numSteps++;
  // } while (t < T - EPS && numSteps <= maxNumSteps);

  // cout << "Total time for computations: " << total_computations / UNIT_SECONDS << " " << LABEL_SECONDS << endl;
  // cout << "Total time for writing: " << total_write / UNIT_SECONDS << " " << LABEL_SECONDS << endl;

  file_sol.close();
  file_E.close();
  return 0;
}