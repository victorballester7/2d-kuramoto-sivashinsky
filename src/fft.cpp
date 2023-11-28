#include "../include/fft.hpp"

#include <cmath>

void fft(double *data, const int n, const int isign, const int normalize) {
  int nn, mmax, m, j, istep, i;
  double wtemp, wr, wpr, wpi, wi, theta, tempr, tempi;
  if (n < 2 || n & (n - 1)) throw("n must be power of 2 in fft");
  nn = n << 1;
  j = 1;
  for (i = 1; i < nn; i += 2) {
    if (j > i) {
      SWAP(data[j - 1], data[i - 1]);
      SWAP(data[j], data[i]);
    }
    m = n;
    while (m >= 2 && j > m) {
      j -= m;
      m >>= 1;
    }
    j += m;
  }
  mmax = 2;
  while (nn > mmax) {
    istep = mmax << 1;
    theta = -isign * (M_PI * 2.0 / mmax);
    wtemp = sin(0.5 * theta);
    wpr = -2.0 * wtemp * wtemp;
    wpi = sin(theta);
    wr = 1.0;
    wi = 0.0;
    for (m = 1; m < mmax; m += 2) {
      for (i = m; i <= nn; i += istep) {
        j = i + mmax;
        tempr = wr * data[j - 1] - wi * data[j];
        tempi = wr * data[j] + wi * data[j - 1];
        data[j - 1] = data[i - 1] - tempr;
        data[j] = data[i] - tempi;
        data[i - 1] += tempr;
        data[i] += tempi;
      }
      wr = (wtemp = wr) * wpr - wi * wpi + wr;
      wi = wi * wpr + wtemp * wpi + wi;
    }
    mmax = istep;
  }
  // From here we have computed the FFT coefficients or the inverse FFT coefficients without any normalization factor. We now apply the normalization factor.

  if (normalize == -1 && isign == -1)
    for (i = 0; i < nn; i++) data[i] /= n;

  if (normalize == 0)
    for (i = 0; i < nn; i++) data[i] /= sqrt(n);

  if (normalize == 1 && isign == 1)
    for (i = 0; i < nn; i++) data[i] /= n;
}

void fft(Vector<double> &data, const int isign, const int normalize) {
  fft(&data[0], data.size() / 2, isign, normalize);
}

void fft(Vector<complex<double>> &data, const int isign, const int normalize) {
  fft((double *)(&data[0]), data.size(), isign, normalize);
}

// void fft_real(Vector<double> &data, const int isign) {
//   int i, i1, i2, i3, i4, n = data.size();
//   double c1 = 0.5, c2, h1r, h1i, h2r, h2i, wr, wi, wpr, wpi, wtemp;
//   // double theta = 3.141592653589793238 / Doub(n >> 1);
//   double theta = M_PI / double(n >> 1);
//   if (isign == 1) {
//     c2 = -0.5;
//     fft(data, 1);
//   } else {
//     c2 = 0.5;
//     theta = -theta;
//   }
//   wtemp = sin(0.5 * theta);
//   wpr = -2.0 * wtemp * wtemp;
//   wpi = sin(theta);
//   wr = 1.0 + wpr;
//   wi = wpi;
//   for (i = 1; i < (n >> 2); i++) {
//     i2 = 1 + (i1 = i + i);
//     i4 = 1 + (i3 = n - i1);
//     h1r = c1 * (data[i1] + data[i3]);
//     h1i = c1 * (data[i2] - data[i4]);
//     h2r = -c2 * (data[i2] + data[i4]);
//     h2i = c2 * (data[i1] - data[i3]);
//     data[i1] = h1r + wr * h2r - wi * h2i;
//     data[i2] = h1i + wr * h2i + wi * h2r;
//     data[i3] = h1r - wr * h2r + wi * h2i;
//     data[i4] = -h1i + wr * h2i + wi * h2r;
//     wr = (wtemp = wr) * wpr - wi * wpi + wr;
//     wi = wi * wpr + wtemp * wpi + wi;
//   }
//   if (isign == 1) {
//     data[0] = (h1r = data[0]) + data[1];
//     data[1] = h1r - data[1];
//   } else {
//     data[0] = c1 * ((h1r = data[0]) + data[1]);
//     data[1] = c1 * (h1r - data[1]);
//     fft(data, -1);
//   }
// }

// void sinfft(Vector<double> &y, const int isign) {
//   int j, n = y.size();
//   double sum, y1, y2, theta, wi = 0.0, wr = 1.0, wpi, wpr, wtemp;
//   theta = M_PI / double(n);
//   wtemp = sin(0.5 * theta);
//   wpr = -2.0 * wtemp * wtemp;
//   wpi = sin(theta);
//   y[0] = 0.0;
//   for (j = 1; j < (n >> 1) + 1; j++) {
//     wr = (wtemp = wr) * wpr - wi * wpi + wr;
//     wi = wi * wpr + wtemp * wpi + wi;
//     y1 = wi * (y[j] + y[n - j]);
//     y2 = 0.5 * (y[j] - y[n - j]);
//     y[j] = y1 + y2;
//     y[n - j] = y1 - y2;
//   }
//   fft_real(y, 1);
//   y[0] *= 0.5;
//   sum = y[1] = 0.0;
//   for (j = 0; j < n - 1; j += 2) {
//     sum += y[j];
//     y[j] = y[j + 1];
//     y[j + 1] = sum;
//   }
//   if (isign == -1) {  // added by Víctor
//     for (int i = 0; i < n; i++) y[i] *= 2.0 / n;
//   }
// }

// void cosft1(Vector<double> &y, const int isign) {
//   int j, n = y.size() - 1;
//   double sum, y1, y2, theta, wi = 0.0, wpi, wpr, wr = 1.0, wtemp;
//   Vector<double> yy(n);
//   theta = M_PI / n;
//   wtemp = sin(0.5 * theta);
//   wpr = -2.0 * wtemp * wtemp;
//   wpi = sin(theta);
//   sum = 0.5 * (y[0] - y[n]);
//   yy[0] = 0.5 * (y[0] + y[n]);
//   for (j = 1; j < n / 2; j++) {
//     wr = (wtemp = wr) * wpr - wi * wpi + wr;
//     wi = wi * wpr + wtemp * wpi + wi;
//     y1 = 0.5 * (y[j] + y[n - j]);
//     y2 = (y[j] - y[n - j]);
//     yy[j] = y1 - wi * y2;
//     yy[n - j] = y1 + wi * y2;
//     sum += wr * y2;
//   }
//   yy[n / 2] = y[n / 2];
//   fft_real(yy, 1);
//   for (j = 0; j < n; j++) y[j] = yy[j];
//   y[n] = y[1];
//   y[1] = sum;
//   for (j = 3; j < n; j += 2) {
//     sum += y[j];
//     y[j] = sum;
//   }
//   if (isign == -1) {  // added by Víctor
//     for (j = 0; j < n + 1; j++) y[j] *= 2.0 / n;
//   }
// }
// void cosft2(Vector<double> &y, const int isign) {
//   int i, n = y.size();
//   double sum, sum1, y1, y2, ytemp, theta, wi = 0.0, wi1, wpi, wpr, wr = 1.0, wr1, wtemp;
//   theta = M_PI_2 / n;
//   wr1 = cos(theta);
//   wi1 = sin(theta);
//   wpr = -2.0 * wi1 * wi1;
//   wpi = sin(2.0 * theta);
//   if (isign == 1) {
//     for (i = 0; i < n / 2; i++) {
//       y1 = 0.5 * (y[i] + y[n - 1 - i]);
//       y2 = wi1 * (y[i] - y[n - 1 - i]);
//       y[i] = y1 + y2;
//       y[n - 1 - i] = y1 - y2;
//       wr1 = (wtemp = wr1) * wpr - wi1 * wpi + wr1;
//       wi1 = wi1 * wpr + wtemp * wpi + wi1;
//     }
//     fft_real(y, 1);
//     for (i = 2; i < n; i += 2) {
//       wr = (wtemp = wr) * wpr - wi * wpi + wr;
//       wi = wi * wpr + wtemp * wpi + wi;
//       y1 = y[i] * wr - y[i + 1] * wi;
//       y2 = y[i + 1] * wr + y[i] * wi;
//       y[i] = y1;
//       y[i + 1] = y2;
//     }
//     sum = 0.5 * y[1];
//     for (i = n - 1; i > 0; i -= 2) {
//       sum1 = sum;
//       sum += y[i];
//       y[i] = sum1;
//     }
//   } else if (isign == -1) {
//     ytemp = y[n - 1];
//     for (i = n - 1; i > 2; i -= 2)
//       y[i] = y[i - 2] - y[i];
//     y[1] = 2.0 * ytemp;
//     for (i = 2; i < n; i += 2) {
//       wr = (wtemp = wr) * wpr - wi * wpi + wr;
//       wi = wi * wpr + wtemp * wpi + wi;
//       y1 = y[i] * wr + y[i + 1] * wi;
//       y2 = y[i + 1] * wr - y[i] * wi;
//       y[i] = y1;
//       y[i + 1] = y2;
//     }
//     fft_real(y, -1);
//     for (i = 0; i < n / 2; i++) {
//       y1 = y[i] + y[n - 1 - i];
//       y2 = (0.5 / wi1) * (y[i] - y[n - 1 - i]);
//       y[i] = 0.5 * (y1 + y2);
//       y[n - 1 - i] = 0.5 * (y1 - y2);
//       wr1 = (wtemp = wr1) * wpr - wi1 * wpi + wr1;
//       wi1 = wi1 * wpr + wtemp * wpi + wi1;
//     }
//   }
// }

void fft_n(double *data, const Vector<int> &nn, const int isign, const int normalize) {
  int idim, i1, i2, i3, i2rev, i3rev, ip1, ip2, ip3, ifp1, ifp2;
  int ibit, k1, k2, n, nprev, nrem, ntot = 1, ndim = nn.size();
  double tempi, tempr, theta, wi, wpi, wpr, wr, wtemp;
  for (idim = 0; idim < ndim; idim++) ntot *= nn[idim];
  if (ntot < 2 || ntot & (ntot - 1)) throw("must have powers of 2 in fourn");
  nprev = 1;
  for (idim = ndim - 1; idim >= 0; idim--) {
    n = nn[idim];
    nrem = ntot / (n * nprev);
    ip1 = nprev << 1;
    ip2 = ip1 * n;
    ip3 = ip2 * nrem;
    i2rev = 0;
    for (i2 = 0; i2 < ip2; i2 += ip1) {
      if (i2 < i2rev) {
        for (i1 = i2; i1 < i2 + ip1 - 1; i1 += 2) {
          for (i3 = i1; i3 < ip3; i3 += ip2) {
            i3rev = i2rev + i3 - i2;
            SWAP(data[i3], data[i3rev]);
            SWAP(data[i3 + 1], data[i3rev + 1]);
          }
        }
      }
      ibit = ip2 >> 1;
      while (ibit >= ip1 && i2rev + 1 > ibit) {
        i2rev -= ibit;
        ibit >>= 1;
      }
      i2rev += ibit;
    }
    ifp1 = ip1;
    while (ifp1 < ip2) {
      ifp2 = ifp1 << 1;
      theta = -isign * M_PI * 2.0 / (ifp2 / ip1);
      wtemp = sin(0.5 * theta);
      wpr = -2.0 * wtemp * wtemp;
      wpi = sin(theta);
      wr = 1.0;
      wi = 0.0;
      for (i3 = 0; i3 < ifp1; i3 += ip1) {
        for (i1 = i3; i1 < i3 + ip1 - 1; i1 += 2) {
          for (i2 = i1; i2 < ip3; i2 += ifp2) {
            k1 = i2;
            k2 = k1 + ifp1;
            tempr = wr * data[k2] - wi * data[k2 + 1];
            tempi = wr * data[k2 + 1] + wi * data[k2];
            data[k2] = data[k1] - tempr;
            data[k2 + 1] = data[k1 + 1] - tempi;
            data[k1] += tempr;
            data[k1 + 1] += tempi;
          }
        }
        wr = (wtemp = wr) * wpr - wi * wpi + wr;
        wi = wi * wpr + wtemp * wpi + wi;
      }
      ifp1 = ifp2;
    }
    nprev *= n;
  }
  // From here we have computed the FFT coefficients or the inverse FFT coefficients without any normalization factor. We now apply the normalization factor.

  if (normalize == -1 && isign == -1)
    for (i1 = 0; i1 < 2 * ntot; i1++) data[i1] /= ntot;

  if (normalize == 0)
    for (i1 = 0; i1 < 2 * ntot; i1++) data[i1] /= sqrt(ntot);

  if (normalize == 1 && isign == 1)
    for (i1 = 0; i1 < 2 * ntot; i1++) data[i1] /= ntot;
}

void fft_n(Vector<double> &data, const Vector<int> &nn, const int isign, const int normalize) {
  fft_n(&data[0], nn, isign, normalize);
}

// void fft_3_real(double *data, double *speq, const int isign,
//                 const int nn1, const int nn2, const int nn3) {
//   int i1, i2, i3, j1, j2, j3, k1, k2, k3, k4;
//   double theta, wi, wpi, wpr, wr, wtemp;
//   double c1, c2, h1r, h1i, h2r, h2i;
//   Vector<int> nn(3);
//   Vector<double *> spq(nn1);
//   for (i1 = 0; i1 < nn1; i1++) spq[i1] = speq + 2 * nn2 * i1;
//   c1 = 0.5;
//   c2 = -0.5 * isign;
//   theta = isign * (M_PI * 2.0 / nn3);
//   wtemp = sin(0.5 * theta);
//   wpr = -2.0 * wtemp * wtemp;
//   wpi = sin(theta);
//   nn[0] = nn1;
//   nn[1] = nn2;
//   nn[2] = nn3 >> 1;
//   if (isign == 1) {
//     fft_n(data, nn, isign, normalize);
//     k1 = 0;
//     for (i1 = 0; i1 < nn1; i1++)
//       for (i2 = 0, j2 = 0; i2 < nn2; i2++, k1 += nn3) {
//         spq[i1][j2++] = data[k1];
//         spq[i1][j2++] = data[k1 + 1];
//       }
//   }
//   for (i1 = 0; i1 < nn1; i1++) {
//     j1 = (i1 != 0 ? nn1 - i1 : 0);
//     wr = 1.0;
//     wi = 0.0;
//     for (i3 = 0; i3 <= (nn3 >> 1); i3 += 2) {
//       k1 = i1 * nn2 * nn3;
//       k3 = j1 * nn2 * nn3;
//       for (i2 = 0; i2 < nn2; i2++, k1 += nn3) {
//         if (i3 == 0) {
//           j2 = (i2 != 0 ? ((nn2 - i2) << 1) : 0);
//           h1r = c1 * (data[k1] + spq[j1][j2]);
//           h1i = c1 * (data[k1 + 1] - spq[j1][j2 + 1]);
//           h2i = c2 * (data[k1] - spq[j1][j2]);
//           h2r = -c2 * (data[k1 + 1] + spq[j1][j2 + 1]);
//           data[k1] = h1r + h2r;
//           data[k1 + 1] = h1i + h2i;
//           spq[j1][j2] = h1r - h2r;
//           spq[j1][j2 + 1] = h2i - h1i;
//         } else {
//           j2 = (i2 != 0 ? nn2 - i2 : 0);
//           j3 = nn3 - i3;
//           k2 = k1 + i3;
//           k4 = k3 + j2 * nn3 + j3;
//           h1r = c1 * (data[k2] + data[k4]);
//           h1i = c1 * (data[k2 + 1] - data[k4 + 1]);
//           h2i = c2 * (data[k2] - data[k4]);
//           h2r = -c2 * (data[k2 + 1] + data[k4 + 1]);
//           data[k2] = h1r + wr * h2r - wi * h2i;
//           data[k2 + 1] = h1i + wr * h2i + wi * h2r;
//           data[k4] = h1r - wr * h2r + wi * h2i;
//           data[k4 + 1] = -h1i + wr * h2i + wi * h2r;
//         }
//       }
//       wr = (wtemp = wr) * wpr - wi * wpi + wr;
//       wi = wi * wpr + wtemp * wpi + wi;
//     }
//   }
//   if (isign == -1) fft_n(data, nn, isign);
// }

// void fft_2_3_real(double *data, double *speq, const int isign,
//                   const int nn1, const int nn2, const int nn3) {
//   int i1, i2, i3, j1, j2, j3, k1, k2, k3, k4;
//   double theta, wi, wpi, wpr, wr, wtemp;
//   double c1, c2, h1r, h1i, h2r, h2i;
//   Vector<int> nn(3);
//   Vector<double *> spq(nn1);
//   for (i1 = 0; i1 < nn1; i1++) spq[i1] = speq + 2 * nn2 * i1;
//   c1 = 0.5;
//   c2 = -0.5 * isign;
//   theta = isign * (M_PI * 2.0 / nn3);
//   wtemp = sin(0.5 * theta);
//   wpr = -2.0 * wtemp * wtemp;
//   wpi = sin(theta);
//   nn[0] = nn1;
//   nn[1] = nn2;
//   nn[2] = nn3 >> 1;
//   if (isign == 1) {
//     fft_n(data, nn, isign);
//     k1 = 0;
//     for (i1 = 0; i1 < nn1; i1++)
//       for (i2 = 0, j2 = 0; i2 < nn2; i2++, k1 += nn3) {
//         spq[i1][j2++] = data[k1];
//         spq[i1][j2++] = data[k1 + 1];
//       }
//   }
//   for (i1 = 0; i1 < nn1; i1++) {
//     j1 = (i1 != 0 ? nn1 - i1 : 0);
//     wr = 1.0;
//     wi = 0.0;
//     for (i3 = 0; i3 <= (nn3 >> 1); i3 += 2) {
//       k1 = i1 * nn2 * nn3;
//       k3 = j1 * nn2 * nn3;
//       for (i2 = 0; i2 < nn2; i2++, k1 += nn3) {
//         if (i3 == 0) {
//           j2 = (i2 != 0 ? ((nn2 - i2) << 1) : 0);
//           h1r = c1 * (data[k1] + spq[j1][j2]);
//           h1i = c1 * (data[k1 + 1] - spq[j1][j2 + 1]);
//           h2i = c2 * (data[k1] - spq[j1][j2]);
//           h2r = -c2 * (data[k1 + 1] + spq[j1][j2 + 1]);
//           data[k1] = h1r + h2r;
//           data[k1 + 1] = h1i + h2i;
//           spq[j1][j2] = h1r - h2r;
//           spq[j1][j2 + 1] = h2i - h1i;
//         } else {
//           j2 = (i2 != 0 ? nn2 - i2 : 0);
//           j3 = nn3 - i3;
//           k2 = k1 + i3;
//           k4 = k3 + j2 * nn3 + j3;
//           h1r = c1 * (data[k2] + data[k4]);
//           h1i = c1 * (data[k2 + 1] - data[k4 + 1]);
//           h2i = c2 * (data[k2] - data[k4]);
//           h2r = -c2 * (data[k2 + 1] + data[k4 + 1]);
//           data[k2] = h1r + wr * h2r - wi * h2i;
//           data[k2 + 1] = h1i + wr * h2i + wi * h2r;
//           data[k4] = h1r - wr * h2r + wi * h2i;
//           data[k4 + 1] = -h1i + wr * h2i + wi * h2r;
//         }
//       }
//       wr = (wtemp = wr) * wpr - wi * wpi + wr;
//       wi = wi * wpr + wtemp * wpi + wi;
//     }
//   }
//   if (isign == -1) fft_n(data, nn, isign);
// }

// void fft_2_3_real(Matrix3d<double> &data, Matrix<double> &speq, const int isign) {
//   if (speq.nrows() != data.dim1() || speq.ncols() != 2 * data.dim2())
//     throw("bad dims in rlft3");
//   fft_3_real(&data[0][0][0], &speq[0][0], isign, data.dim1(), data.dim2(), data.dim3());
// }

// void fft_2_3_real(Matrix<double> &data, Vector<double> &speq, const int isign) {
//   if (speq.size() != unsigned(2) * data.nrows()) throw("bad dims in rlft3");
//   fft_3_real(&data[0][0], &speq[0], isign, 1, data.nrows(), data.ncols());
// }
