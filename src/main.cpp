// -------------------------------------------------
// -------------------------------------------------
// -------------------- HEADER ---------------------
// -------------------------------------------------
// -------------------------------------------------

#include <fftw3.h>

#include <chrono>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "../include/energy.hpp"
#include "../include/misc.hpp"

using namespace std;

#define START_TIMER() begin = chrono::steady_clock::now()
#define END_TIMER() end = chrono::steady_clock::now()
#define ADD_TIME_TO(x) x += chrono::duration_cast<chrono::microseconds>(end - begin).count()
#define FFTW_FLAG FFTW_ESTIMATE  // options: FFTW_ESTIMATE, FFTW_MEASURE.
// FFTW_MEASURE instructs FFTW to run and measure the execution time of several FFTs in order to find the best way to compute the transform of size n.This process takes some time(usually a few seconds), depending on your machine and on the size of the transform.FFTW_ESTIMATE, on the contrary, does not run any computation and just builds a reasonable plan that is probably sub - optimal.In short, if your program performs many transforms of the same size and initialization time is not important, use FFTW_MEASURE; otherwise use the estimate.

#define WRITE_SOL() ((write_solution) && (numSteps % freq_anim_sol == 0))
#define WRITE_FREQ() ((write_freq) && (numSteps % freq_anim_freq == 0))
#define WRITE_E() ((write_energy) && (numSteps % freq_plot_e == 0))

#define IMEXEULER 1
#define IMEXBDF2 2
#define IMEXRK4 3

#define NONLINEAR_COMPUTATIONS()                                           \
  {                                                                        \
    /* Step in finite differences */                                       \
    /* Multiply Fu by i * kx and store the result in aux_x_fourier */      \
    /* Multiply Fu by i * ky and store the result in aux_y_fourier */      \
    for (uint i = 0; i < dim_f; i++) {                                     \
      aux_x_fourier[i][0] = -Fu[i][1] * kx[i];                             \
      aux_x_fourier[i][1] = Fu[i][0] * kx[i];                              \
      aux_y_fourier[i][0] = -Fu[i][1] * ky[i];                             \
      aux_y_fourier[i][1] = Fu[i][0] * ky[i];                              \
    }                                                                      \
    /* Perform inverse discrete Fourier transform (DFT) */                 \
    fftw_execute(p_back_aux_x);                                            \
    fftw_execute(p_back_aux_y);                                            \
    /* Normalize the inverse DFT and Compute (u_x)^2 and (u_y)^2 */        \
    /* We divide by dim because fftw will not normalize the next FFT */    \
    for (uint i = 0; i < dim; i++) {                                       \
      aux_x[i] = (aux_x[i]) * (aux_x[i]) / dim;                            \
      aux_y[i] = (aux_y[i]) * (aux_y[i]) / dim;                            \
    }                                                                      \
    /* Compute the DFT of (u_x)^2 and store the result in aux_x_fourier */ \
    fftw_execute(p_for_aux_x);                                             \
    /* Compute the DFT of (u_y)^2 and store the result in aux_y_fourier */ \
    fftw_execute(p_for_aux_y);                                             \
  }

#define NONLINEAR_COMPUTATIONS_2()                                         \
  {                                                                        \
    /* Step in finite differences */                                       \
    /* Multiply Fu by i * kx and store the result in aux_x_fourier */      \
    /* Multiply Fu by i * ky and store the result in aux_y_fourier */      \
    for (uint i = 0; i < dim_f; i++) {                                     \
      aux_x_fourier_2[i][0] = -Fu_2[i][1] * kx[i];                         \
      aux_x_fourier_2[i][1] = Fu_2[i][0] * kx[i];                          \
      aux_y_fourier_2[i][0] = -Fu_2[i][1] * ky[i];                         \
      aux_y_fourier_2[i][1] = Fu_2[i][0] * ky[i];                          \
    }                                                                      \
    /* Perform inverse discrete Fourier transform (DFT) */                 \
    fftw_execute(p_back_aux_x_2);                                          \
    fftw_execute(p_back_aux_y_2);                                          \
    /* Normalize the inverse DFT and Compute (u_x)^2 and (u_y)^2 */        \
    for (uint i = 0; i < dim; i++) {                                       \
      aux_x[i] = (aux_x[i]) * (aux_x[i]) / dim;                            \
      aux_y[i] = (aux_y[i]) * (aux_y[i]) / dim;                            \
    }                                                                      \
    /* Compute the DFT of (u_x)^2 and store the result in aux_x_fourier */ \
    fftw_execute(p_for_aux_x_2);                                           \
    /* Compute the DFT of (u_y)^2 and store the result in aux_y_fourier */ \
    fftw_execute(p_for_aux_y_2);                                           \
  }

#define NONLINEAR_TERM(i, c) (-0.5 * (aux_x_fourier[(i)][(c)] + nu2 / nu1 * aux_y_fourier[(i)][(c)]))
#define NONLINEAR_TERM_2(i, c) (-0.5 * (aux_x_fourier_2[(i)][(c)] + nu2 / nu1 * aux_y_fourier_2[(i)][(c)]))
#define INTEGRAL_NONLINEAR_TERM (0.5 * (aux_x_fourier[0][0] + nu2 / nu1 * aux_y_fourier[0][0]))
#define INTEGRAL_NONLINEAR_TERM_2 (0.5 * (aux_x_fourier_2[0][0] + nu2 / nu1 * aux_y_fourier_2[0][0]))

#define STEP_IMEXEULER()                                                        \
  {                                                                             \
    /* Create the nonlinear term */                                             \
    NONLINEAR_COMPUTATIONS();                                                   \
    /* Update the iterate of the finite differences scheme */                   \
    for (uint i = 0; i < dim_f; i++) {                                          \
      Fu[i][0] = C1[i] * ((1 + dt * c) * Fu[i][0] + dt * NONLINEAR_TERM(i, 0)); \
      Fu[i][1] = C1[i] * ((1 + dt * c) * Fu[i][1] + dt * NONLINEAR_TERM(i, 1)); \
    }                                                                           \
    /* Increase the time */                                                     \
    t += dt;                                                                    \
  }

#define STEP_IMEXEULER_NONHOMO()                                                                \
  {                                                                                             \
    /* Create the nonlinear term */                                                             \
    NONLINEAR_COMPUTATIONS();                                                                   \
    /* Non-homogeneous term */                                                                  \
    Fw[0][0] = 1 + cos(t);                                                                      \
    /* Compute the DFT of the non-homogeneous term */                                           \
    /* Update the iterate of the finite differences scheme */                                   \
    for (uint i = 0; i < dim_f; i++) {                                                          \
      Fu[i][0] = C1[i] * ((1 + dt * c) * Fu[i][0] + dt * NONLINEAR_TERM(i, 0) + dt * Fw[i][0]); \
      Fu[i][1] = C1[i] * ((1 + dt * c) * Fu[i][1] + dt * NONLINEAR_TERM(i, 1) + dt * Fw[i][1]); \
    }                                                                                           \
    /* Increase the time */                                                                     \
    t += dt;                                                                                    \
  }

#define STEP_IMEXBDF2()                                                                                                                                 \
  {                                                                                                                                                     \
    /* Fu_2 contains the n-th iterate of the finite differences scheme */                                                                               \
    /* Fu contains the (n+1)-th iterate of the finite differences scheme */                                                                             \
    /* We want to compute the (n+2)-th iterate of the finite differences scheme */                                                                      \
                                                                                                                                                        \
    /* Create the nonlinear term */                                                                                                                     \
    NONLINEAR_COMPUTATIONS();                                                                                                                           \
    NONLINEAR_COMPUTATIONS_2();                                                                                                                         \
    /* Update the iterate of the finite differences scheme */                                                                                           \
    for (uint i = 0; i < dim_f; i++) {                                                                                                                  \
      Fu_aux[i][0] = C2[i] * (Fu[i][0] * 2 * (1 + dt * c) - Fu_2[i][0] * (0.5 + dt * c) + 2 * dt * NONLINEAR_TERM(i, 0) - dt * NONLINEAR_TERM_2(i, 0)); \
      Fu_aux[i][1] = C2[i] * (Fu[i][1] * 2 * (1 + dt * c) - Fu_2[i][1] * (0.5 + dt * c) + 2 * dt * NONLINEAR_TERM(i, 1) - dt * NONLINEAR_TERM_2(i, 1)); \
    }                                                                                                                                                   \
    memcpy(Fu_2, Fu, dim_f * sizeof(fftw_complex));                                                                                                     \
    memcpy(Fu, Fu_aux, dim_f * sizeof(fftw_complex));                                                                                                   \
    /* Increase the time */                                                                                                                             \
    t += dt;                                                                                                                                            \
  }

#define STEP_IMEXBDF2_NONHOMO()                                                                                                                                                    \
  {                                                                                                                                                                                \
    /* Fu_2 contains the n-th iterate of the finite differences scheme */                                                                                                          \
    /* Fu contains the (n+1)-th iterate of the finite differences scheme */                                                                                                        \
    /* We want to compute the (n+2)-th iterate of the finite differences scheme */                                                                                                 \
                                                                                                                                                                                   \
    /* Create the nonlinear term */                                                                                                                                                \
    NONLINEAR_COMPUTATIONS();                                                                                                                                                      \
    NONLINEAR_COMPUTATIONS_2();                                                                                                                                                    \
    /* Update the iterate of the finite differences scheme */                                                                                                                      \
    Fw[0][0] = 1 + cos(t - dt);                                                                                                                                                    \
    Fww[0][0] = 1 + cos(t);                                                                                                                                                        \
    for (uint i = 0; i < dim_f; i++) {                                                                                                                                             \
      Fu_aux[i][0] = C2[i] * (Fu[i][0] * 2 * (1 + dt * c) - Fu_2[i][0] * (0.5 + dt * c) + 2 * dt * (NONLINEAR_TERM(i, 0) + Fww[i][0]) - dt * (NONLINEAR_TERM_2(i, 0) + Fw[i][0])); \
      Fu_aux[i][1] = C2[i] * (Fu[i][1] * 2 * (1 + dt * c) - Fu_2[i][1] * (0.5 + dt * c) + 2 * dt * (NONLINEAR_TERM(i, 1) + Fww[i][1]) - dt * (NONLINEAR_TERM_2(i, 1) + Fw[i][1])); \
    }                                                                                                                                                                              \
    memcpy(Fu_2, Fu, dim_f * sizeof(fftw_complex));                                                                                                                                \
    memcpy(Fu, Fu_aux, dim_f * sizeof(fftw_complex)); /* Increase the time */                                                                                                      \
    t += dt;                                                                                                                                                                       \
  }
#define STEP_IMEXRK4()                                                                                                                                              \
  {                                                                                                                                                                 \
    /*computation of the 1st stage*/                                                                                                                                \
    NONLINEAR_COMPUTATIONS();                                                                                                                                       \
    for (uint i = 0; i < dim_f; i++) {                                                                                                                              \
      Fu_2[i][0] = 1. / denom[i * 4] * (Fu[i][0] + dt * (betaI[0] * L[i] * Fu[i][0] + alphaE[0] * NONLINEAR_TERM(i, 0)));                                           \
      Fu_2[i][1] = 1. / denom[i * 4] * (Fu[i][1] + dt * (betaI[0] * L[i] * Fu[i][1] + alphaE[0] * NONLINEAR_TERM(i, 1)));                                           \
    }                                                                                                                                                               \
                                                                                                                                                                    \
    /*computation of the 2nd stage*/                                                                                                                                \
    NONLINEAR_COMPUTATIONS_2();                                                                                                                                     \
    for (uint i = 0; i < dim_f; i++) {                                                                                                                              \
      Fu[i][0] = 1. / denom[i * 4 + 1] * (Fu_2[i][0] + dt * (betaI[1] * L[i] * Fu_2[i][0] + alphaE[1] * NONLINEAR_TERM_2(i, 0) + betaE[1] * NONLINEAR_TERM(i, 0))); \
      Fu[i][1] = 1. / denom[i * 4 + 1] * (Fu_2[i][1] + dt * (betaI[1] * L[i] * Fu_2[i][1] + alphaE[1] * NONLINEAR_TERM_2(i, 1) + betaE[1] * NONLINEAR_TERM(i, 1))); \
    }                                                                                                                                                               \
                                                                                                                                                                    \
    /*computation of the 3rd stage*/                                                                                                                                \
    NONLINEAR_COMPUTATIONS();                                                                                                                                       \
    for (uint i = 0; i < dim_f; i++) {                                                                                                                              \
      Fu_2[i][0] = 1. / denom[i * 4 + 2] * (Fu[i][0] + dt * (betaI[2] * L[i] * Fu[i][0] + alphaE[2] * NONLINEAR_TERM(i, 0) + betaE[2] * NONLINEAR_TERM_2(i, 0)));   \
      Fu_2[i][1] = 1. / denom[i * 4 + 2] * (Fu[i][1] + dt * (betaI[2] * L[i] * Fu[i][1] + alphaE[2] * NONLINEAR_TERM(i, 1) + betaE[2] * NONLINEAR_TERM_2(i, 1)));   \
    }                                                                                                                                                               \
                                                                                                                                                                    \
    /*computation of the 4th stage*/                                                                                                                                \
    NONLINEAR_COMPUTATIONS_2();                                                                                                                                     \
    /* here Fu = Fu_2*/                                                                                                                                             \
    for (uint i = 0; i < dim_f; i++) {                                                                                                                              \
      Fu[i][0] = 1. / denom[i * 4 + 3] * (Fu_2[i][0] + dt * (betaI[3] * L[i] * Fu_2[i][0] + alphaE[3] * NONLINEAR_TERM_2(i, 0) + betaE[3] * NONLINEAR_TERM(i, 0))); \
      Fu[i][1] = 1. / denom[i * 4 + 3] * (Fu_2[i][1] + dt * (betaI[3] * L[i] * Fu_2[i][1] + alphaE[3] * NONLINEAR_TERM_2(i, 1) + betaE[3] * NONLINEAR_TERM(i, 1))); \
    }                                                                                                                                                               \
    /* Increase the time */                                                                                                                                         \
    t += dt;                                                                                                                                                        \
  }

#define STEP_IMEXRK4_NONHOMO()                                                                                                                                                                          \
  {                                                                                                                                                                                                     \
    /*computation of the 1st stage*/                                                                                                                                                                    \
    Fw[0][0] = 1 + cos(t);                                                                                                                                                                              \
    NONLINEAR_COMPUTATIONS();                                                                                                                                                                           \
    for (uint i = 0; i < dim_f; i++) {                                                                                                                                                                  \
      Fu_2[i][0] = 1. / denom[i * 4] * (Fu[i][0] + dt * (betaI[0] * L[i] * Fu[i][0] + alphaE[0] * NONLINEAR_TERM(i, 0) + alphaE[0] * Fw[i][0]));                                                        \
      Fu_2[i][1] = 1. / denom[i * 4] * (Fu[i][1] + dt * (betaI[0] * L[i] * Fu[i][1] + alphaE[0] * NONLINEAR_TERM(i, 1) + alphaE[0] * Fw[i][1]));                                                        \
    }                                                                                                                                                                                                   \
                                                                                                                                                                                                        \
    /*computation of the 2nd stage*/                                                                                                                                                                    \
    Fw[0][0] = 1 + cos(t + dt * 14.0 / 25.0);                                                                                                                                                           \
    NONLINEAR_COMPUTATIONS_2();                                                                                                                                                                         \
    for (uint i = 0; i < dim_f; i++) {                                                                                                                                                                  \
      Fu[i][0] = 1. / denom[i * 4 + 1] * (Fu_2[i][0] + dt * (betaI[1] * L[i] * Fu_2[i][0] + alphaE[1] * NONLINEAR_TERM_2(i, 0) + betaE[1] * NONLINEAR_TERM(i, 0) + (betaE[1] + alphaE[1]) * Fw[i][0])); \
      Fu[i][1] = 1. / denom[i * 4 + 1] * (Fu_2[i][1] + dt * (betaI[1] * L[i] * Fu_2[i][1] + alphaE[1] * NONLINEAR_TERM_2(i, 1) + betaE[1] * NONLINEAR_TERM(i, 1) + (betaE[1] + alphaE[1]) * Fw[i][1])); \
    }                                                                                                                                                                                                   \
                                                                                                                                                                                                        \
    /*computation of the 3rd stage*/                                                                                                                                                                    \
    Fw[0][0] = 1 + cos(t + dt * 41.0 / 50.0);                                                                                                                                                           \
    NONLINEAR_COMPUTATIONS();                                                                                                                                                                           \
    for (uint i = 0; i < dim_f; i++) {                                                                                                                                                                  \
      Fu_2[i][0] = 1. / denom[i * 4 + 2] * (Fu[i][0] + dt * (betaI[2] * L[i] * Fu[i][0] + alphaE[2] * NONLINEAR_TERM(i, 0) + betaE[2] * NONLINEAR_TERM_2(i, 0) + (betaE[2] + alphaE[2]) * Fw[i][0]));   \
      Fu_2[i][1] = 1. / denom[i * 4 + 2] * (Fu[i][1] + dt * (betaI[2] * L[i] * Fu[i][1] + alphaE[2] * NONLINEAR_TERM(i, 1) + betaE[2] * NONLINEAR_TERM_2(i, 1) + (betaE[2] + alphaE[2]) * Fw[i][1]));   \
    }                                                                                                                                                                                                   \
                                                                                                                                                                                                        \
    /*computation of the 4th stage*/                                                                                                                                                                    \
    NONLINEAR_COMPUTATIONS_2();                                                                                                                                                                         \
    /* here Fu = Fu_2*/                                                                                                                                                                                 \
    for (uint i = 0; i < dim_f; i++) {                                                                                                                                                                  \
      Fu[i][0] = 1. / denom[i * 4 + 3] * (Fu_2[i][0] + dt * (betaI[3] * L[i] * Fu_2[i][0] + alphaE[3] * NONLINEAR_TERM_2(i, 0) + betaE[3] * NONLINEAR_TERM(i, 0) + (betaE[3] + alphaE[3]) * Fw[i][0])); \
      Fu[i][1] = 1. / denom[i * 4 + 3] * (Fu_2[i][1] + dt * (betaI[3] * L[i] * Fu_2[i][1] + alphaE[3] * NONLINEAR_TERM_2(i, 1) + betaE[3] * NONLINEAR_TERM(i, 1) + (betaE[3] + alphaE[3]) * Fw[i][1])); \
    }                                                                                                                                                                                                   \
    /* Increase the time */                                                                                                                                                                             \
    t += dt;                                                                                                                                                                                            \
  }

#define CLEANUP()                         \
  {                                       \
    if (write_solution) file_sol.close(); \
    if (write_energy) file_E.close();     \
    /* Clean variables */                 \
    free(kx);                             \
    free(ky);                             \
    free(C1);                             \
    free(C2);                             \
    free(denom);                          \
    free(L);                              \
    /* FFT Cleanup */                     \
    fftw_destroy_plan(p_for_u);           \
    fftw_destroy_plan(p_back_u);          \
    fftw_destroy_plan(p_for_aux_x);       \
    fftw_destroy_plan(p_for_aux_y);       \
    fftw_destroy_plan(p_back_aux_x);      \
    fftw_destroy_plan(p_back_aux_y);      \
    fftw_destroy_plan(p_for_aux_x_2);     \
    fftw_destroy_plan(p_for_aux_y_2);     \
    fftw_destroy_plan(p_back_aux_x_2);    \
    fftw_destroy_plan(p_back_aux_y_2);    \
    fftw_free(u);                         \
    fftw_free(aux_x);                     \
    fftw_free(aux_y);                     \
    fftw_free(Fu);                        \
    fftw_free(Fu_2);                      \
    fftw_free(Fu_aux);                    \
    fftw_free(aux_x_fourier);             \
    fftw_free(aux_y_fourier);             \
    fftw_free(aux_x_fourier_2);           \
    fftw_free(aux_y_fourier_2);           \
    fftw_cleanup();                       \
  }

#define EPS 1e-8

const double alphaI[4] = {343038331393.0 / 1130875731271.0, 288176579239.0 / 1140253497719.0, 253330171251.0 / 677500478386.0, 189462239225.0 / 1091147436423.0};
const double betaI[4] = {35965327958.0 / 140127563663.0, 19632212512.0 / 2700543775099.0, -173747147147.0 / 351772688865.0, 91958533623.0 / 727726057489.0};
const double alphaE[4] = {14.0 / 25.0, 777974228744.0 / 1346157007247.0, 251277807242.0 / 1103637129625.0, 113091689455.0 / 220187950967.0};
const double betaE[4] = {0.0, -251352885992.0 / 790610919619.0, -383714262797.0 / 1103637129625.0, -403360439203.0 / 1888264787188.0};

// -------------------------------------------------
// -------------------------------------------------
// ----------------- MAIN FUNCTION -----------------
// -------------------------------------------------
// -------------------------------------------------

int main(void) {
  // ----------------- PARAMETERS ------------------
  const string filename_solution = "data/solution.txt";            // name of the output file to write the solution
  const string filename_energy = "data/energy.txt";                // name of the output file to write the energy
  const string filename_energy_return = "data/energy_return.txt";  // name of the output file to write the return map of the energy
  const string filename_input = "data/input.txt";                  // name of the input file
  const string filename_freq = "data/freq.txt";                    // name of the output file to write the frequency of the energy
  const string space = "    ";                                     // space to print
  const uint per = 10;                                             // progress percentage interval to print (each count%)
  uint count = per;
  chrono::steady_clock::time_point begin, end;                                                                                                      // variables to measure the time
  int64_t total_write = 0, total_fftw_computations = 0, total_fft_precomputation = 0, total_energy_computations = 0, total_other_computations = 0;  // variables to measure the time
  int nx, ny;                                                                                                                                       // number of points in x and y (must be a power of 2), do NOT use uint fftw uses int
  double nu1, nu2;                                                                                                                                  // parameters of the system (by default)
  double dt;                                                                                                                                        // initial step size
  double T;                                                                                                                                         // final time of integration
  double cutoff_time;                                                                                                                               // time at which it starts the stationary regime
  // DO NOT DECREASE THE NU'S BELOW 0.2 BECAUSE THE SOLUTION WILL START TO EXPLODE
  bool averaged_solution;  // 1 if we want to average the solution, 0 otherwise
  bool write_solution;     // whether to write the solution to a file or not
  bool write_energy;       // whether to write the energy to a file or not
  bool write_freq;         // whether to write the frequency of the energy to a file or not
  uint freq_anim_sol;      // frequency of the plot of the solution
  uint freq_plot_e;        // frequency of the plot of the energy
  uint freq_anim_freq;     // frequency of the plot of the frequency of the energy
  uint method;
  string method_name;
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
    file_input >> tmp >> dt;
    file_input >> tmp >> T;
    file_input >> tmp >> cutoff_time;
    file_input >> tmp >> write_solution;
    file_input >> tmp >> write_freq;
    file_input >> tmp >> write_energy;
    file_input >> tmp >> averaged_solution;
    file_input >> tmp >> freq_anim_sol;
    file_input >> tmp >> freq_anim_freq;
    file_input >> tmp >> freq_plot_e;
    file_input >> tmp >> method;
  }
  file_input.close();
  // Alternative for debugging
  // nx = 32;
  // ny = 32;
  // nu1 = 1.0;
  // nu2 = 1.0;
  // dt = 5.0e-3;
  // T = 5.0;
  // cutoff_time = 100.0;
  // write_solution = true;
  // write_energy = true;
  // averaged_solution = true;
  // freq_plot_sol = 1;
  // freq_plot_e = 1;
  // method = IMEXEULER;
  if (method == IMEXRK4) {
    method_name = "IMEX-RK4";
  } else if (method == IMEXBDF2) {
    method_name = "IMEX-BDF2";
  } else if (method == IMEXEULER) {
    method_name = "IMEX-Euler";
  } else {
    cout << "Method not recognized" << endl;
    return 1;
  }
  // -----------------------------------------------

  // plot the parameters
  cout << "Size of the grid:               " << space << nx << " x " << ny << endl;
  cout << "Parameters:                     " << space << "nu1 = " << nu1 << endl;
  cout << "                                " << space << "nu2 = " << nu2 << endl;
  cout << "Final time:                     " << space << T << endl;
  cout << "Step size:                      " << space << dt << endl;
  cout << "Use averaged solution?          " << space << (averaged_solution ? "yes" : "no") << endl;
  cout << "Write physical solution to file?" << space << (write_solution ? "yes" : "no") << endl;
  cout << "Write Fourier solution to file? " << space << (write_freq ? "yes" : "no") << endl;
  cout << "Write energy to file?           " << space << (write_energy ? "yes" : "no") << endl;
  cout << "Method:                         " << space << method_name << endl;
  END_TIMER();
  ADD_TIME_TO(total_write);

  // ----------------- PARAMETERS ------------------
  int ny_complex = ny / 2 + 1;           // number of points in y in the frequency space to optimize the memory in the FFT
  uint dim = (uint)(nx * ny);            // dimension of the physical system
  uint dim_f = (uint)(nx * ny_complex);  // dimension of the system in the frequency space
  double fraction_completed = T / 100.;  // fraction of the integration time to print
  double t = 0.;                         // initial time of integration
  uint numSteps = 0;                     // number of steps
  // vector to store all the solutions
  // double energy = 0, energy2 = 0, energy3 = 0, dE = 0, dE_1 = 0, dE_2 = 0, tn, En;  // energy of the system
  My_double energy = {0, 0}, energy2 = {0, 0}, energy3 = {0, 0};  // energy of the system
  My_double dE = {0, 0}, dE_1 = {0, 0}, dE_2 = {0, 0};            // derivative of the energy
  double factor_E = (2 * M_PI) * (2 * M_PI) / (nx * ny);          // factor to compute the energy
  double factor_dE = (2 * M_PI) * (2 * M_PI) / (dt * nx * ny);    // factor to compute the derivative of the energy
  double tn, En;
  bool En_changed = false;
  double c = 1.0 + 1.0 / (nu1);  // constant term to ensure stability of the scheme (see: p227 in Nonlinear dynamics of surfactant–laden multilayer shear flows and related systems)
  double *kx = (double *)malloc(dim_f * sizeof(double));
  double *ky = (double *)malloc(dim_f * sizeof(double));
  double *C1 = (double *)malloc(dim_f * sizeof(double));
  double *C2 = (double *)malloc(dim_f * sizeof(double));
  // -----------------------------------------------

  // ----------------- FFT SETUP -------------------
  double *u = fftw_alloc_real(dim);   // u is the solution,
  double *v = fftw_alloc_real(dim);   // v is the solution,
  double *w = fftw_alloc_real(dim);   // w is the solution,
  double *ww = fftw_alloc_real(dim);  // w is the solution,
  double *aux_x = fftw_alloc_real(dim);
  double *aux_y = fftw_alloc_real(dim);
  fftw_complex *Fu = fftw_alloc_complex(dim_f);   // Fu is the DFT of u
  fftw_complex *Fv = fftw_alloc_complex(dim_f);   // Fv is the DFT of v
  fftw_complex *Fw = fftw_alloc_complex(dim_f);   // Fw is the DFT of w
  fftw_complex *Fww = fftw_alloc_complex(dim_f);  // FU is the DFT of u
  fftw_complex *Fu_2 = fftw_alloc_complex(dim_f);
  fftw_complex *Fu_aux = fftw_alloc_complex(dim_f);
  fftw_complex *aux_x_fourier = fftw_alloc_complex(dim_f);
  fftw_complex *aux_y_fourier = fftw_alloc_complex(dim_f);
  fftw_complex *aux_x_fourier_2 = fftw_alloc_complex(dim_f);
  fftw_complex *aux_y_fourier_2 = fftw_alloc_complex(dim_f);
  // ----------------------------------------

  fftw_plan p_for_u, p_back_u, p_for_aux_x, p_for_aux_y, p_back_aux_x, p_back_aux_y, p_back_aux_x_2, p_back_aux_y_2, p_for_aux_x_2, p_for_aux_y_2, p_for_w, p_for_ww;

  START_TIMER();
  // Remember that for inverse transforms, that is for _c2r_ transforms, the input array is overwritten, so you should copy it before if you want to preserve the input data.
  // Also remember that the fftw computes unnormalized transforms, so doing a forward transform followed by a backward transform will multiply the input by n.
  p_for_u = fftw_plan_dft_r2c_2d(nx, ny, u, Fu, FFTW_FLAG);       // forward
  p_back_u = fftw_plan_dft_c2r_2d(nx, ny, Fu_aux, u, FFTW_FLAG);  // backward

  // auxiliar
  p_for_w = fftw_plan_dft_r2c_2d(nx, ny, w, Fw, FFTW_FLAG);     // forward
  p_for_ww = fftw_plan_dft_r2c_2d(nx, ny, ww, Fww, FFTW_FLAG);  // forward

  p_for_aux_x = fftw_plan_dft_r2c_2d(nx, ny, aux_x, aux_x_fourier, FFTW_FLAG);   // forward
  p_for_aux_y = fftw_plan_dft_r2c_2d(nx, ny, aux_y, aux_y_fourier, FFTW_FLAG);   // forward
  p_back_aux_x = fftw_plan_dft_c2r_2d(nx, ny, aux_x_fourier, aux_x, FFTW_FLAG);  // backward
  p_back_aux_y = fftw_plan_dft_c2r_2d(nx, ny, aux_y_fourier, aux_y, FFTW_FLAG);  // backward

  p_for_aux_x_2 = fftw_plan_dft_r2c_2d(nx, ny, aux_x, aux_x_fourier_2, FFTW_FLAG);   // forward
  p_for_aux_y_2 = fftw_plan_dft_r2c_2d(nx, ny, aux_y, aux_y_fourier_2, FFTW_FLAG);   // forward
  p_back_aux_x_2 = fftw_plan_dft_c2r_2d(nx, ny, aux_x_fourier_2, aux_x, FFTW_FLAG);  // backward
  p_back_aux_y_2 = fftw_plan_dft_c2r_2d(nx, ny, aux_y_fourier_2, aux_y, FFTW_FLAG);  // backward

  END_TIMER();
  ADD_TIME_TO(total_fft_precomputation);
  // -----------------------------------------------

  START_TIMER();
  set_data(u, nx, ny);
  set_wave_numbers(kx, ky, nx, ny_complex);
  set_C1(C1, kx, ky, nx, ny_complex, dt, nu1, nu2, c);
  set_C2(C2, kx, ky, nx, ny_complex, dt, nu1, nu2, c);
  fftw_execute(p_for_u);  // Now Fu contains the DFT of u
  for (uint i = 0; i < dim_f; i++) {
    Fu[i][0] /= dim;
    Fu[i][1] /= dim;
  }
  if (method == IMEXBDF2) {
    memcpy(Fu_2, Fu, dim_f * sizeof(fftw_complex));
  }
  // ---------------- SCHEME SETUP -----------------
  // only needed for method == IMEXRK4
  double *denom = (double *)malloc(dim_f * 4 * sizeof(double));
  double *L = (double *)malloc(dim_f * sizeof(double));
  double aux;
  if (method == IMEXRK4) {
    for (uint i = 0; i < dim_f; i++) {
      aux = kx[i] * kx[i] + (nu2 / nu1) * ky[i] * ky[i];
      L[i] = aux * (1.0 - nu1 * aux);
      for (uint j = 0; j < 4; j++) denom[i * 4 + j] = 1.0 - dt * L[i] * alphaI[j];
    }
  }
  // -----------------------------------------------
  END_TIMER();
  ADD_TIME_TO(total_other_computations);

  // file handling
  START_TIMER();
  ofstream file_sol, file_E, file_En, file_freq;  // output files
  // if (averaged_solution && (WRITE_SOL() || WRITE_E())) {  // Fu[0], the fourier coefficient with k_1 = k_2 = 0, is the mean of the solution
  //   for (uint i = 0; i < dim; i++) u[i] -= Fu[0][0];
  // }
  if (WRITE_SOL()) {
    file_sol.open(filename_solution);
    write(u, nx, ny, t, file_sol);  // write the solution into the file
  }
  if (WRITE_E()) {  // write the energy into the file
    file_E.open(filename_energy);
    file_En.open(filename_energy_return);
    // in the 4th column plot u(pi, pi)
    file_E << t << " " << (energy.int_part + energy.frac_part) << " " << (dE.int_part + dE.frac_part) << " " << u[(ny / 2) * (nx + 1)] << endl;  // write the energy into the file
  }
  if (WRITE_FREQ()) {
    file_freq.open(filename_freq);
    // write_bis(Fu, nx, ny_complex, t, file_freq);
    write_fourier(Fu, nx, ny_complex, t, file_freq);
  }
  END_TIMER();
  ADD_TIME_TO(total_write);
  numSteps++;
  // fftw_execute(p_for_aux_x_2);
  // for (uint i = 0; i < nx; i++) {
  //   for (uint j = 0; j < ny; j++) {
  //     w[i * ny + j] = g(i * 2 * M_PI / nx, j * 2 * M_PI / ny, t) / dim;
  //     ww[i * ny + j] = g(i * 2 * M_PI / nx, j * 2 * M_PI / ny, t) / dim;
  //   }
  // }                       /* Compute the DFT of the non-homogeneous term */
  // fftw_execute(p_for_w);  /* Update the iterate of the finite differences scheme */
  // fftw_execute(p_for_ww); /* Update the iterate of the finite differences scheme */
  // dt = 0.1;
  // int numDt = 12;
  // double anterior = 0, aux_aux;
  // memcpy(Fv, Fu, dim_f * sizeof(fftw_complex));
  // memcpy(v, u, dim * sizeof(double));
  // for (int k = 0; k < numDt; k++) {
  //   memcpy(u, v, dim * sizeof(double));
  //   memcpy(Fu, Fv, dim_f * sizeof(fftw_complex));
  //   memcpy(Fu_2, Fu, dim_f * sizeof(fftw_complex));
  //   // restart other parameters:
  //   set_C1(C1, kx, ky, nx, ny_complex, dt, nu1, nu2, c);
  //   set_C2(C2, kx, ky, nx, ny_complex, dt, nu1, nu2, c);
  //   t = 0;
  //   if (method == IMEXRK4) {
  //     for (uint i = 0; i < dim_f; i++) {
  //       aux_aux = kx[i] * kx[i] + (nu2 / nu1) * ky[i] * ky[i];
  //       L[i] = aux_aux * (1.0 - nu1 * aux_aux);
  //       for (uint j = 0; j < 4; j++) denom[i * 4 + j] = 1.0 - dt * L[i] * alphaI[j];
  //     }
  //   }
  //   numSteps = 1;
  //   anterior = aux;
  double max;
  do {
    memcpy(v, u, dim * sizeof(double));
    START_TIMER();
    if (method == IMEXEULER) {
      STEP_IMEXEULER();
    } else if (method == IMEXBDF2) {
      if (numSteps == 1) {
        STEP_IMEXEULER();
      } else {
        STEP_IMEXBDF2();
      }
    } else if (method == IMEXRK4) {
      STEP_IMEXRK4();
    }
    // normalize the inverse DFT. We first substract the mean of the solution and then we normalize it
    // we would lose the information of Fu, so we save it in Fu_aux and now we will lose the information of Fu_aux
    memcpy(Fu_aux, Fu, dim_f * sizeof(fftw_complex));
    fftw_execute(p_back_u);  // Now u contains the inverse DFT of Fu, which is the new iterate of the solution in the physical space
                             // if (averaged_solution) {  // Fu[0], the fourier coefficient with k_1 = k_2 = 0, is the mean of the solution
    if (averaged_solution)
      for (uint i = 0; i < dim; i++) u[i] -= Fu[0][0];  // we need to normalize Fu accordingly in order to be the mean of the solution
    // else
    //   for (uint i = 0; i < dim; i++) u[i] /= dim;  // we need to normalize Fu accordingly in order to be the mean of the solution
    // }
    END_TIMER();
    ADD_TIME_TO(total_fftw_computations);
    // compute the energy (we have to do it each time because if not we would get strange derivatives of the energy, we would be a wrong 'energy2' when writing the derivative)
    START_TIMER();
    energy3 = energy2;
    energy2 = energy;
    energy = My_E(u, nx, ny);
    END_TIMER();
    ADD_TIME_TO(total_energy_computations);
    // energy = sqrt(energy);
    START_TIMER();
    dE_2 = dE_1;
    dE_1 = dE;
    dE.int_part = (energy.int_part - energy2.int_part);
    dE.frac_part = (energy.frac_part - energy2.frac_part);
    aux = (dE.int_part + dE.frac_part) * (dE_1.int_part + dE_1.frac_part);
    if (aux < 0 && aux < -EPS) {  // then we have a zero in the derivative between t and t - dt
      tn = lagrange_root(t, dt, dE_2, dE_1, dE);
      En = lagrange_eval(tn, t, dt, energy3, energy2, energy) * factor_E;
      En_changed = true;
    } else
      En_changed = false;
    END_TIMER();

    // for (uint i = 0; i < dim_f; i++) {
    //   aux_x_fourier[i][0] = -Fu[i][1] * kx[i];
    //   aux_x_fourier[i][1] = Fu[i][0] * kx[i];
    //   aux_y_fourier[i][0] = -Fu[i][1] * ky[i];
    //   aux_y_fourier[i][1] = Fu[i][0] * ky[i];
    // } /* Perform inverse discrete Fourier transform (DFT) */

    // fftw_execute(p_back_aux_x);
    // fftw_execute(p_back_aux_y);

    // // file_sol << t << " " << max(aux_x, dim) << " " << max(aux_y, dim) << endl;
    // for (uint i = 0; i < dim; i++) {
    //   aux_x[i] = aux_x[i] * aux_x[i] + aux_y[i] * aux_y[i];
    // }
    max = 0;
    for (uint i = 0; i < dim; i++) {
      aux_x[i] = abs((u[i] - v[i]) / u[i]);
      if (aux_x[i] > max) max = aux_x[i];
    }

    ADD_TIME_TO(total_other_computations);

    START_TIMER();
    if (WRITE_SOL()) {
      write(u, nx, ny, t, file_sol);  // write the solution into the file
      // write(aux_x, nx, ny, t, file_sol);  // write the solution into the file
      // file_sol << t << "    " << max << endl;
    }
    // we plot the energy En at each time independently of freq_plot_e
    if (write_energy && En_changed) file_En << tn << " " << En << endl;  // write the energy into the file

    if (WRITE_E()) {  // write the energy into the file
      // plot more digits of the energy at each time
      file_E << t << " " << (energy.int_part + energy.frac_part) * factor_E << " " << (dE.int_part + dE.frac_part) * factor_dE << " " << u[(ny / 2) * (nx + 1)] << endl;  // energy, derivative of the energy, u(pi, pi)
    }

    if (WRITE_FREQ()) {
      // write_bis(Fu, nx, ny_complex, t, file_freq);
      write_fourier(Fu, nx, ny_complex, t, file_freq);
    }

    if (t > fraction_completed * count - EPS) {
      cout << count << "%" << endl;
      count += per;
    }
    END_TIMER();
    ADD_TIME_TO(total_write);
    numSteps++;
  } while (t < T - EPS);  // We add 1 because we count the initial condition as a step
  // compute the fft of U
  // memcpy(Fu_aux, Fu, dim_f * sizeof(fftw_complex));
  // fftw_execute(p_back_u);
  // for (uint i = 0; i < nx; i++) {
  //   for (uint j = 0; j < ny; j++) {
  //     u[i * ny + j] -= U(i * 2 * M_PI / nx, j * 2 * M_PI / ny, t);
  //   }
  // }
  // aux = l2(u, nx, ny);
  // cout << "l2 norm of the error (dt = " << dt << "): " << aux;
  //   if (k > 0) {
  //     if (method == IMEXEULER)
  //       cout << " ratio = " << aux / anterior << " log ratio = " << log(anterior / aux) / log(2) << endl;
  //     else if (method == IMEXBDF2)
  //       cout << " ratio = " << aux / (anterior) << " log ratio = " << log(anterior / aux) / log(2) << endl;
  //     else if (method == IMEXRK4)
  //       cout << " ratio = " << aux / (anterior) << " log ratio = " << log(anterior / aux) / log(2) << endl;
  //   } else {
  //     cout << " ratio = 0" << endl;
  //   }
  //   dt /= 2.;
  // }
  // fftw_execute(p_for_aux_x_2);

  // normalize fft
  // for (uint i = 0; i < dim_f; i++) {
  //   aux_x_fourier_2[i][0] /= dim;
  //   aux_x_fourier_2[i][1] /= dim;
  // }

  // creation of tmp files to comunicate with bash scripts in order to plot or not the solution and the energy
  START_TIMER();
  ofstream file_tmp;
  if (write_solution || write_freq) {
    file_tmp.open("data/tmp_write_sol.txt");
    if (write_solution && !write_freq)
      file_tmp << 1 << endl;  // only plot of u
    else if (write_freq && !write_solution)
      file_tmp << 2 << endl;  // only plot of Fu
    else
      file_tmp << 3 << endl;  // plot both
    file_tmp.close();
  }
  if (write_energy) {
    file_tmp.open("data/tmp_write_E.txt");
    file_tmp << cutoff_time << endl;  // we do this to pass the cutoff time to the bash scriptand then to the python script
    file_tmp.close();
  }
  END_TIMER();
  ADD_TIME_TO(total_write);

  print("Total time for fftw precomputation: ", total_fft_precomputation);
  print("Total time for fftw computations: ", total_fftw_computations);
  print("Total time for energy computations: ", total_energy_computations);
  print("Total time for other computations: ", total_other_computations);
  print("Total time for read/write: ", total_write);

  CLEANUP();

  return 0;
}
