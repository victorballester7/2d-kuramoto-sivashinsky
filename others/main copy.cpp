// #include "../include/field.h"
// #include "../include/rk78.h"

#include "../include/LinearAlgebra.hpp"

// x should contain all the initial conditions (x0, y0, x1, y1, x2, y2, ...) of the different orbits
// the result will be stored in result which will be a vector of the form
// [x0(t0), y0(t0), x0(t1), y0(t1), ..., x0(tn), y0(tn), x1(t0), y1(t0), x1(t1), y1(t1), ..., x1(tn), y1(tn), ...]

// double* integration(int dim, int numSteps, int num_orbits, double x[], double h, double hmin, double hmax, double tol_rk78, int maxNumStepsFlow, double mu, double xmin, double xmax, double ymin, double ymax) {
//   double* result = malloc(sizeof(double) * (numSteps + 1) * dim * num_orbits);
//   double t = 0;
//   double T = hmax;
//   Args arg = {malloc(sizeof(double)), xmin, xmax, ymin, ymax};
//   arg.mu[0] = mu;
//   // print x
//   // for (int i = 0; i < num_orbits; i++) {
//   //   printf("x = %lf, y = %lf\n", x[i * dim], x[i * dim + 1]);
//   // }
//   for (int i = 0; i < num_orbits; i++) {
//     t = 0;
//     memcpy(result + i * dim * (numSteps + 1), x + i * dim, sizeof(double) * dim);  // we copy the initial conditions
//     // printf("Hola-i=%i\n", i);
//     for (int j = 1; j < numSteps + 1; j++) {
//       // printf("Hola-j1=%i\n", j);
//       // printf("hola1, t = %lf, h = %lf, hmin = %lf, hmax = %lf\n", t, h, hmin, hmax);
//       // print x
//       // if (i == 1) {
//       //   printf("x[%i] = %lf, y[%i] = %lf\n", i, x[i * dim], i, x[i * dim + 1]);
//       //   printf("mu=%lf, t=%lf, x[%i] = %lf, y[%i] = %lf, h = %lf, T=%lf, hmin = %lf, hmax = %lf, tol = %lf, n = %i,maxNumSteps=%d\n", arg.mu[0], t, i, x[i * dim], i, x[i * dim + 1], h, T, hmin, hmax, tol_rk78, dim, maxNumStepsFlow);
//       // }
//       // printf("x[%i] = %lf, y[%i] = %lf\n", i, x[i * dim], i, x[i * dim + 1]);
//       // printf("result[%i] = %lf, result[%i] = %lf\n", (j - 1) * dim + i * dim * (numSteps + 1), result[(j - 1) * dim + i * dim * (numSteps + 1)], (j - 1) * dim + i * dim * (numSteps + 1) + 1, result[(j - 1) * dim + i * dim * (numSteps + 1) + 1]);
//       if (flow(&t, x + i * dim, &h, T, hmin, hmax, tol_rk78, maxNumStepsFlow, dim, field, &arg)) {
//         printf("Error in integration\n");
//         // printf("Hola-error-j=%i\n", j);
//         exit(1);
//       }
//       // printf("x[%i] = %lf, y[%i] = %lf\n", i, x[i * dim], i, x[i * dim + 1]);
//       // printf("result[%i] = %lf, result[%i] = %lf\n", (j - 1) * dim + i * dim * (numSteps + 1), result[(j - 1) * dim + i * dim * (numSteps + 1)], (j - 1) * dim + i * dim * (numSteps + 1) + 1, result[(j - 1) * dim + i * dim * (numSteps + 1) + 1]);
//       // compare if the new iterate is the same as the previous one
//       if ((x[i * dim] - result[(j - 1) * dim + i * dim * (numSteps + 1)] < tol_rk78 && x[i * dim] - result[(j - 1) * dim + i * dim * (numSteps + 1)] > -tol_rk78) && (x[i * dim + 1] - result[(j - 1) * dim + i * dim * (numSteps + 1) + 1] < tol_rk78 && x[i * dim + 1] - result[(j - 1) * dim + i * dim * (numSteps + 1) + 1] > -tol_rk78)) {
//         // printf("Error new iterate is the same as the previous one-i= %i, j = %i\n", i, j);
//         for (int k = j; k < numSteps + 1; k++) {  // we copy nan value to the rest of the vector to save time
//           result[k * dim + i * dim * (numSteps + 1)] = NAN;
//           result[k * dim + i * dim * (numSteps + 1) + 1] = NAN;
//         }
//         break;
//       }
//       // printf("hola2, t = %lf, h = %lf, hmin = %lf, hmax = %lf\n", t, h, hmin, hmax);
//       memcpy(result + j * dim + i * dim * (numSteps + 1), x + i * dim, sizeof(double) * dim);
//       // printf("Hola-j2=%i\n", j);
//       // for (int k = 0; k < dim; k++) {
//       //   printf("result[%i] = %lf     ", k, result[j * dim + i * dim * numSteps + k]);
//       // }
//       // printf("\n");
//     }
//   }

//   // for (int i = 0; i < num_orbits; i++) {
//   //   for (int j = 0; j < numSteps + 1; j++) {
//   //     printf("x=%lf, y=%lf\n", result[j * dim + i * dim * numSteps], result[j * dim + i * dim * numSteps + 1]);
//   //     printf("\n");
//   //   }
//   //   printf("\n---------------------------------------\n");
//   // }
//   return result;
// }

int main(void) {
  Vector<double> v(3);
  v[0] = 1;
  v[1] = 2;
  v[2] = 3;

  v.print();

  // double* x = malloc(sizeof(double) * 2);
  // x[0] = 0.5;
  // x[1] = -0.25;
  // int n = 2;
  // double t = 0;
  // double T = -0.5;
  // double h = -0.02;
  // double hmin = -0.01;
  // double hmax = -0.04;
  // double tol = 1e-8;
  // int numMax = 100;
  // int numSteps = 1000;
  // Args arg = {malloc(sizeof(double)), -2, 2, -2, 2};
  // arg.mu[0] = -0.1;
  // double* result = integration(n, numSteps, 1, x, h, hmin, hmax, tol, numMax, arg.mu[0], arg.xmin, arg.xmax, arg.ymin, arg.ymax);
  // // double* result = integration(n, numSteps, 1, x, h, hmin, hmax, tol, numMax, arg.mu[0], arg.xmin, arg.xmax, arg.ymin, arg.ymax);
  // // print result for copy to np.array
  // printf("L = np.array([\n");
  // for (int i = 0; i < 100; i++) {
  //   printf("[%lf, %lf],\n", result[i * 2], result[i * 2 + 1]);
  // }
  // printf("])\n");
  // for (int i = 0; i < numSteps; i++) {
  //   // printf("Hola-i=%i\n", i);
  //   printf("mu= %lf, t=%lf, x[0] = %lf, x[1] = %lf, h = %lf, T=%lf, hmin = %lf, hmax = %lf, tol = %lf, n = %i,maxNumSteps=%d\n", arg.mu[0], t, x[0], x[1], h, T, hmin, hmax, tol, n, numMax);
  //   flow(&t, x, &h, T, hmin, hmax, tol, numMax, n, field, &arg);
  //   for (int k = 0; k < n; k++) {
  //     printf("x[%i] = %lf     ", k, x[k]);
  //   }
  //   printf("\n");
  // }
  return 0;
}

// compilate with (the include files are in include/):
// gcc -Wall -Iinclude/ src/integration.c src/rk78.c src/field.c -o bin/integration -lm
