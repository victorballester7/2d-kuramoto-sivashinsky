#ifndef FIELD_HPP
#define FIELD_HPP

#include "LinearAlgebra.hpp"
#include "fft.hpp"

typedef struct args {
  double nu1, nu2;
  Vector<int> nn;
  Vector<double> k1, k2;
  Vector<double> tmp;  // tmp = (k₁² +(ν₂ / ν₁) * k₂²)
} Args;

// -----------------------------------------------------
// field
// -----------------------------------------------------
// Purpose:
// 	Computes the vector field
// -----------------------------------------------------
int field(int n, double t, double x[], double f[], void *param) {
  // ∂ₜû = -1 / 2 * ℱ([ℱ⁻¹(i * k₁ * û)]² + (ν₂ / ν₁) * [ℱ⁻¹(i * k₂ * û)]²)
  //      +(k₁² +(ν₂ / ν₁) * k₂²) * û
  //      - ν₁ *(k₁² +(ν₂ / ν₁) * k₂²)² * û
  Args *args = (Args *)param;

  Vector<double> X(x, n);

  Vector<double> tmp1(n), tmp2(n), aux(n);

  // try with throwing an exception
  try {
    tmp1 = args->k1 * X;
    tmp2 = args->k2 * X;

    for (int i = 0; i < n; i += 2) {  // multiply by i
      tmp1[i] = -tmp1[i + 1];
      tmp1[i + 1] = tmp1[i];
      tmp2[i] = -tmp2[i + 1];
      tmp2[i + 1] = tmp2[i];
    }

    fft_n(tmp1, args->nn, -1, 1);
    fft_n(tmp2, args->nn, -1, 1);

    aux = tmp1 * tmp1 + (args->nu2 / args->nu1) * tmp2 * tmp2;

    fft_n(aux, args->nn, 1, 1);

    aux *= -0.5;
    aux += args->tmp * (1.0 - args->nu1 * args->tmp) * X;

    // copy aux to f
    for (int i = 0; i < n; i++) f[i] = aux[i];
  } catch (const std::exception &e) {
    std::cout << "Problem in evaluating the function: " << e.what() << std::endl;
    return 1;
  }

  return 0;
}

#endif  // FIELD_HPP
