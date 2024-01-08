#include "../include/scheme_fd.hpp"

// void stepFiniteDiff(vector<double> &x, double &t, const Args &prm) {
//   // 2D Kuramoto - Sivashinsky equation
//   //  ∂ₜu + 1 / 2 * | ∇u |² + Δu + Δ²u = 0
//   // Doing a change of variable to rescale the domain to [0, 2π] we get the equation:
//   //  ∂ₜu + 1 / 2 * ((∂ₓu)² + (ν₂ / ν₁) * (∂ᵧu)²) + [∂ₓ²u + (ν₂ / ν₁) * ∂ᵧ²u] + ν₁ *[∂ₓ⁴u + 2 * (ν₂ / ν₁) * ∂ₓ²∂ᵧ²u + (ν₂ / ν₁)² *∂ᵧ⁴u] = 0
//   // Taking the Fourier transform of the equation we get:
//   //  ∂ₜû = -1 / 2 * ℱ([ℱ⁻¹(i * k₁ * û)]² + (ν₂ / ν₁) * [ℱ⁻¹(i * k₂ * û)]²)
//   //        +(k₁² +(ν₂ / ν₁) * k₂²) * û
//   //        - ν₁ *(k₁² +(ν₂ / ν₁) * k₂²)² * û
//   // We can write this as:
//   //  ∂ₜû = C * û + f(û)
//   // where C = -(k₁² +(ν₂ / ν₁) * k₂²) + ν₁ *(k₁² +(ν₂ / ν₁) * k₂²)² is the linear part
//   // and f(û) = -1 / 2 * ℱ([ℱ⁻¹(i * k₁ * û)]² + (ν₂ / ν₁) * [ℱ⁻¹(i * k₂ * û)]²) is the nonlinear part

//   // We can solve this equation using the linear-implicit-nonlinear-explicit Euler method
//   // û_{n+1} - û_n = h * (C * û_{n+1} + f(û_n))
//   // û_{n+1} = (I - h * C)^-1 * (û_n + h * f(û_n))

//   // Implicit-explicit Euler method
//   // vector<double> aux1, aux2;
//   // double tmp;
//   // aux1 = prm.k1 * x;
//   // aux2 = prm.k2 * x;
//   // for (int i = 0; i < prm.k1.size(); i += 2) {  // multiply by i
//   //   tmp = aux1[i];
//   //   aux1[i] = -aux1[i + 1];
//   //   aux1[i + 1] = tmp;
//   //   tmp = aux2[i];
//   //   aux2[i] = -aux2[i + 1];
//   //   aux2[i + 1] = tmp;
//   // }
//   // fft_n(aux1, prm.nn, -1, 1);
//   // fft_n(aux2, prm.nn, -1, 1);
//   // aux1 *= aux1;
//   // aux2 *= aux2;
//   // fft_n(aux1, prm.nn, 1, 1);
//   // fft_n(aux2, prm.nn, 1, 1);

//   x = prm.tmp * (x + (prm.h * N(x, prm)));
//   // x = prm.tmp * x;

//   // only linear part
//   // x = prm.tmp * x;
//   t += prm.h;
// }

// // vector<double> N(vector<double> &x, const Args &prm) {
// //   // nonlinear part of the equation
// //   vector<double> aux1, aux2;
// //   double tmp;
// //   aux1 = prm.k1 * x;
// //   aux2 = prm.k2 * x;
// //   for (int i = 0; i < prm.k1.size(); i += 2) {  // multiply by i
// //     tmp = aux1[i];
// //     aux1[i] = -aux1[i + 1];
// //     aux1[i + 1] = tmp;
// //     tmp = aux2[i];
// //     aux2[i] = -aux2[i + 1];
// //     aux2[i + 1] = tmp;
// //   }
// //   fft_n(aux1, prm.nn, -1, 1);
// //   fft_n(aux2, prm.nn, -1, 1);
// //   aux1 *= aux1;
// //   aux2 *= aux2;
// //   fft_n(aux1, prm.nn, 1, 1);
// //   fft_n(aux2, prm.nn, 1, 1);

// //   return -0.5 * (aux1 + ((prm.nu2 / prm.nu1) * aux2));
// // }

// void stepIMEXRK4(vector<double> &x, double &t, const Args &prm) {
//   // IMEXRK4 scheme for the 2D Kuramoto - Sivashinsky equation
//   // suppose we have the equation:
//   //  ∂ₜu = L * u + N(u)
//   // where L is the linear part and N(u) is the nonlinear part
//   // The 4-steps schemes writes as:
//   //  û⁽¹⁾ = ûⁿ + h * (ɑ_1I * L * û⁽¹⁾ + β_1I * L * ûⁿ + ɑ_1E * N(ûⁿ) + 0)
//   //  û⁽²⁾ = û⁽¹⁾ + h * (ɑ_2I * L * û⁽²⁾ + β_2I * L * û⁽¹⁾ + ɑ_2E * N(û⁽¹⁾) + β_2E * N(ûⁿ))
//   //  û⁽³⁾ = û⁽²⁾ + h * (ɑ_3I * L * û⁽³⁾ + β_3I * L * û⁽²⁾ + ɑ_3E * N(û⁽²⁾) + β_3E * N(û⁽¹⁾))
//   //  û⁽⁴⁾ = û⁽³⁾ + h * (ɑ_4I * L * û⁽⁴⁾ + β_4I * L * û⁽³⁾ + ɑ_4E * N(û⁽³⁾) + β_4E * N(û⁽²⁾))

//   // Hence we have the following equations for the intermediate steps:
//   //  û⁽¹⁾ = (1 - h * ɑ_1I * L)^-1 * (ûⁿ + h * (β_1I * L * ûⁿ + ɑ_1E * N(ûⁿ)))
//   //  û⁽²⁾ = (1 - h * ɑ_2I * L)^-1 * (û⁽¹⁾ + h * (β_2I * L * û⁽¹⁾ + ɑ_2E * N(û⁽¹⁾) + β_2E * N(ûⁿ)))
//   //  û⁽³⁾ = (1 - h * ɑ_3I * L)^-1 * (û⁽²⁾ + h * (β_3I * L * û⁽²⁾ + ɑ_3E * N(û⁽²⁾) + β_3E * N(û⁽¹⁾)))
//   //  û⁽⁴⁾ = (1 - h * ɑ_4I * L)^-1 * (û⁽³⁾ + h * (β_4I * L * û⁽³⁾ + ɑ_4E * N(û⁽³⁾) + β_4E * N(û⁽²⁾)))

//   vector<double> u1, u2, u3, u4;
//   vector<double> denom[4];
//   for (int i = 0; i < 4; i++) {
//     denom[i] = 1. - ((prm.h * alphaI[i]) * prm.L);
//   }
//   u1 = 1. / denom[0] * (x + prm.h * ((betaI[0] * prm.L * x) + (alphaE[0] * N(x, prm))));
//   u2 = 1. / denom[1] * (u1 + prm.h * ((betaI[1] * prm.L * u1) + (alphaE[1] * N(u1, prm)) + (betaE[1] * N(x, prm))));
//   u3 = 1. / denom[2] * (u2 + prm.h * ((betaI[2] * prm.L * u2) + (alphaE[2] * N(u2, prm)) + (betaE[2] * N(u1, prm))));
//   u4 = 1. / denom[3] * (u3 + prm.h * ((betaI[3] * prm.L * u3) + (alphaE[3] * N(u3, prm)) + (betaE[3] * N(u2, prm))));
//   x = u4;
//   t += prm.h;
// }