#ifndef FIELD_HPP
#define FIELD_HPP

#include "LinearAlgebra.hpp"
#include "fft.hpp"

typedef struct args {
  double nu1, nu2;
  vector<uint> nn;
  Vector<double> k1, k2;
  Vector<double> tmp;  // tmp = (k₁² +(ν₂ / ν₁) * k₂²)
} Args;

Vector<double> g(Vector<double> &x) {
  Vector<double> y(x.size());
  for (int i = 0; i < x.size(); i++) {
    y[i] = x[i] * x[i] * x[i] * x[i];
  }
  return y;
}

// -----------------------------------------------------
// field
// -----------------------------------------------------
// Purpose:
// 	Computes the vector field
// -----------------------------------------------------
int field(int n, double t, double x[], double f[], void *param) {
  // 2D Kuramoto - Sivashinsky equation
  //  ∂ₜu + 1 / 2 * | ∇u |² + Δu + Δ²u = 0
  // Doing a change of variable to rescale the domain to[0, 2π] we get the equation:
  //  ∂ₜu + 1 / 2 * ((∂ₓu)² + (ν₂ / ν₁) * (∂ᵧu)²) + [∂ₓ²u + (ν₂ / ν₁) * ∂ᵧ²u] + ν₁ *[∂ₓ⁴u + 2 * (ν₂ / ν₁) * ∂ₓ²∂ᵧ²u + (ν₂ / ν₁)² *∂ᵧ⁴u] = 0
  // Taking the Fourier transform of the equation we get:
  //  ∂ₜû = -1 / 2 * ℱ([ℱ⁻¹(i * k₁ * û)]² + (ν₂ / ν₁) * [ℱ⁻¹(i * k₂ * û)]²)
  //        +(k₁² +(ν₂ / ν₁) * k₂²) * û
  //        - ν₁ *(k₁² +(ν₂ / ν₁) * k₂²)² * û
  Args *args = (Args *)param;

  Vector<double> X(x, n);

  Vector<double> tmp1(n), tmp2(n), aux(n);

  // try with throwing an exception
  try {
    // tmp1 = args->k1 * X;
    // tmp2 = args->k2 * X;

    // for (int i = 0; i < n; i += 2) {  // multiply by i
    //   tmp1[i] = -tmp1[i + 1];
    //   tmp1[i + 1] = tmp1[i];
    //   tmp2[i] = -tmp2[i + 1];
    //   tmp2[i + 1] = tmp2[i];
    // }

    // fft_n(tmp1, args->nn, -1, 1);
    // fft_n(tmp2, args->nn, -1, 1);

    // aux = tmp1 * tmp1 + (args->nu2 / args->nu1) * tmp2 * tmp2;

    // fft_n(aux, args->nn, 1, 1);

    // aux *= -0.5;
    // aux += args->tmp * (1.0 - args->nu1 * args->tmp) * X;

    // // copy aux to f
    // for (int i = 0; i < n; i++) f[i] = aux[i];

    // we first try with the heat equation: ∂ₜu = Δu + g(u)
    // Fourier transform: ∂ₜû = -(k₁² + k₂²) * û + ℱ(g(ℱ⁻¹(û)))
    tmp1 = -1.0 * (args->k1 * args->k1 + args->k2 * args->k2) * X;

    // aux = X;
    // fft_n(aux, args->nn, -1, 1);

    // aux = g(aux);

    // fft_n(aux, args->nn, 1, 1);

    // tmp1 += aux;

    // copy tmp1 to f
    for (int i = 0; i < n; i++) f[i] = tmp1[i];

  } catch (const std::exception &e) {
    std::cout << "Problem in evaluating the function: " << e.what() << std::endl;
    return 1;
  }

  return 0;
}

#endif  // FIELD_HPP
