#include "../include/setup.hpp"

#include <fstream>

#include "../include/fft.hpp"

double u0(double x, double y) {  // initial condition
  if (x < 0) x += 2 * M_PI;      // periodic extension
  if (y < 0) y += 2 * M_PI;      // periodic extension

  // return sin(x) + sin(2 * y) + cos(7 * exp(-x));
  return sin(x) + sin(y) + sin(x + y);
}

// double heat_kernel(double t, double x, double y) {
//   return 1.0 / (4 * M_PI * t) * exp(-(x * x + y * y) / (4 * t));
// }

// double u(double t, double x, double y, int nx, int ny) {
//   // integrate numerically f * H, where H is the heat kernel
//   double xmin = 0, xmax = 2 * M_PI;
//   double ymin = 0, ymax = 2 * M_PI;
//   double hx = (xmax - xmin) / nx;
//   double hy = (ymax - ymin) / ny;
//   double sum = 0.0;

//   for (double i = xmin; i < xmax; i += hx) {
//     for (double j = ymin; j < ymax; j += hy) {
//       sum += f(i, j) * heat_kernel(t, x - i, y - j);
//     }
//   }
//   return sum * hx * hy;
// }

double write(Vector<double> &x, Vector<int> nn, double t, ofstream &file, bool print_complex) {  // write into file
  file << t << endl;
#define round_to_zero(x) (fabs(x) < 1e-10 ? 0 : x)

  for (int i = 0; i < nn[0]; i++) {
    for (int j = 0; j < 2 * nn[1] - 2; j += 2) {
      if (print_complex)
        file << round_to_zero(x[2 * i * nn[1] + j]) << showpos << round_to_zero(x[2 * i * nn[1] + j + 1]) << noshowpos << "i ";
      else
        file << round_to_zero(x[2 * i * nn[1] + j]) << " ";
    }
    if (print_complex)
      file << round_to_zero(x[2 * i * nn[1] + 2 * nn[1] - 2]) << showpos << round_to_zero(x[2 * i * nn[1] + 2 * nn[1] - 1]) << noshowpos << "i" << endl;
    else
      file << round_to_zero(x[2 * i * nn[1] + 2 * nn[1] - 2]) << endl;
  }
  file << endl;
#undef round_to_zero
  return 0;
}

void stepFiniteDiff(Vector<double> &x, double &t, double h, Args &prm) {
  // 2D Kuramoto - Sivashinsky equation
  //  ∂ₜu + 1 / 2 * | ∇u |² + Δu + Δ²u = 0
  // Doing a change of variable to rescale the domain to[0, 2π] we get the equation:
  //  ∂ₜu + 1 / 2 * ((∂ₓu)² + (ν₂ / ν₁) * (∂ᵧu)²) + [∂ₓ²u + (ν₂ / ν₁) * ∂ᵧ²u] + ν₁ *[∂ₓ⁴u + 2 * (ν₂ / ν₁) * ∂ₓ²∂ᵧ²u + (ν₂ / ν₁)² *∂ᵧ⁴u] = 0
  // Taking the Fourier transform of the equation we get:
  //  ∂ₜû = -1 / 2 * ℱ([ℱ⁻¹(i * k₁ * û)]² + (ν₂ / ν₁) * [ℱ⁻¹(i * k₂ * û)]²)
  //        +(k₁² +(ν₂ / ν₁) * k₂²) * û
  //        - ν₁ *(k₁² +(ν₂ / ν₁) * k₂²)² * û
  // We can write this as:
  //  ∂ₜû = C * û + f(û)
  // where C = -(k₁² +(ν₂ / ν₁) * k₂²) + ν₁ *(k₁² +(ν₂ / ν₁) * k₂²)² is the linear part
  // and f(û) = -1 / 2 * ℱ([ℱ⁻¹(i * k₁ * û)]² + (ν₂ / ν₁) * [ℱ⁻¹(i * k₂ * û)]²) is the nonlinear part

  // We can solve this equation using the linear-implicit-nonlinear-explicit Euler method
  // û_{n+1} - û_n = h * (C * û_{n+1} + f(û_n))
  // û_{n+1} = (I + h * C)^-1 * (û_n + h * f(û_n))

  Vector<double> aux1 = x, aux2 = x;
  double tmp;
  aux1 = prm.k1 * x;
  aux2 = prm.k2 * x;
  for (int i = 0; i < prm.k1.size(); i += 2) {  // multiply by i
    tmp = aux1[i];
    aux1[i] = -aux1[i + 1];
    aux1[i + 1] = tmp;
    tmp = aux2[i];
    aux2[i] = -aux2[i + 1];
    aux2[i + 1] = tmp;
  }
  fft_n(aux1, prm.nn, -1, 1);
  fft_n(aux2, prm.nn, -1, 1);
  aux1 *= aux1;
  aux2 *= aux2;
  fft_n(aux1, prm.nn, 1, 1);
  fft_n(aux2, prm.nn, 1, 1);

  x = prm.tmp * (x - h * 0.5 * (aux1 + ((prm.nu2 / prm.nu1) * aux2)));
  t += h;
}

void set_data(Vector<double> &x, Args &prm) {
  double a1 = 2 * M_PI / prm.nx;
  double a2 = 2 * M_PI / prm.ny;
  for (int i = 0; i < prm.nx; i++) {
    for (int j = 0; j < 2 * prm.ny; j += 2) {
      x[2 * i * prm.ny + j] = u0(prm.k1[2 * i * prm.ny + j] * a1, prm.k2[2 * i * prm.ny + j] * a2);
      x[2 * i * prm.ny + j + 1] = 0;
    }
  }
}