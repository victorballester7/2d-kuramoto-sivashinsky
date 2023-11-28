#ifndef FFT_HPP
#define FFT_HPP

#include <bits/stdc++.h>

#include "LinearAlgebra.hpp"
using namespace std;

template <class T>
inline void SWAP(T &a, T &b) {
  T dum = a;
  a = b;
  b = dum;
}

/// @brief Compute the FFT (if isign = 1) or the IFFT (if isign = -1) of a complex vector data of length 2 * n.
/// @param data Vector of data to be transformed, sorted as (real0, imag0, real1, imag1, ...). It will be replaced by the result.
/// @param n Length of the vector data and must be a power of 2.
/// @param isign If isign is 1, compute the FFT, if isign is -1, compute the IFFT.
/// @param normalize The are three commonly used definitions of the Fourier transform:
//    If normalize = 0, the factor 1 / sqrt(n) is included in the forward and inverse transform.
//    If normalize = 1, the factor 1 / n is included in the forward transform (that is, the fourier coefficients get normalized).
//    If normalize = -1, the factor 1 / n is included in the inverse transform (that is, the inverse fourier coefficients get normalized).
// Note: for consistent results, the same value of normalize should be used for both the forward and inverse transform using the same data.
void fft(double *data, const int n, const int isign, const int normalize = 1);

/// @brief Compute the FFT (if isign = 1) or the IFFT (if isign = -1) of a real vector data of length 2 * n.
/// @param data Vector of data to be transformed, sorted as (real0, real1, ...).
/// @param n Length of the vector data and must be a power of 2.
/// @param isign If isign is 1, compute the FFT, if isign is -1, compute the IFFT.
/// @param normalize The are three commonly used definitions of the Fourier transform:
//    If normalize = 0, the factor 1 / sqrt(n) is included in the forward and inverse transform.
//    If normalize = 1, the factor 1 / n is included in the forward transform (that is, the fourier coefficients get normalized).
//    If normalize = -1, the factor 1 / n is included in the inverse transform (that is, the inverse fourier coefficients get normalized).
// Note: for consistent results, the same value of normalize should be used for both the forward and inverse transform using the same data.
void fft(Vector<double> &data, const int isign, const int normalize = 1);

/// @brief Compute the FFT (if isign = 1) or the IFFT (if isign = -1) of a complex vector data of length n.
/// @param data Vector of data to be transformed.
/// @param isign If isign is 1, compute the FFT, if isign is -1, compute the IFFT.
/// @param normalize The are three commonly used definitions of the Fourier transform:
//    If normalize = 0, the factor 1 / sqrt(n) is included in the forward and inverse transform.
//    If normalize = 1, the factor 1 / n is included in the forward transform (that is, the fourier coefficients get normalized).
//    If normalize = -1, the factor 1 / n is included in the inverse transform (that is, the inverse fourier coefficients get normalized).
// Note: for consistent results, the same value of normalize should be used for both the forward and inverse transform using the same data.
void fft(Vector<complex<double>> &data, const int isign, const int normalize = 1);

// // void fft_real(Vector<double> &data, const int isign);

// /// @brief Compute the sine fast transform of a real vector data of length n.
// /// @param data Vector of data to be transformed. It will be replaced by the result. The length of data must be a power of 2.
// /// @param isign If isign is 1, compute the sine fast transform, if isign is -1, compute the inverse sine fast transform.
// void sinft(Vector<double> &y, const int isign);

// /// @brief Compute the cosine fast transform of a real vector data of length n + 1, where n must be a power of 2.
// /// @param data Vector of data to be transformed. It will be replaced by the result.
// /// @param isign If isign is 1, compute the cosine fast transform, if isign is -1, compute the inverse cosine fast transform.
// void cost1(Vector<double> &y, const int isign);

// void cost2(Vector<double> &y, const int isign);

/// @brief Compute the ndim-dimensional FFT (if isign = 1) or the ndim-dimensional IFFT (if isign = -1) of a complex vector data.
/// @param data Vector of data to be transformed. It has the length of nn[0] * nn[1] * nn[2] * ...  and will be replaced by the result at the end of the function. The complex data is stored as (real0, imag0, real1, imag1, ...). For a two-dimensional array, this is equivalent to storing the array by rows.
/// @param nn Vector of integers specifying the lengths of each dimension. All nn[i] must be powers of 2.
/// @param isign If isign is 1, compute the FFT, if isign is -1, compute the IFFT.
/// @param normalize The are three commonly used definitions of the Fourier transform:
//    If normalize = 0, the factor 1 / sqrt(n) is included in the forward and inverse transform.
//    If normalize = 1, the factor 1 / n is included in the forward transform (that is, the fourier coefficients get normalized).
//    If normalize = -1, the factor 1 / n is included in the inverse transform (that is, the inverse fourier coefficients get normalized).
// Note: for consistent results, the same value of normalize should be used for both the forward and inverse transform using the same data.
// The normalize factor used in numpy would be equivalent to set normalize = -1.
void fft_n(double *data, const Vector<int> &nn, const int isign, const int normalize = 1);

/// @brief Compute the ndim-dimensional FFT (if isign = 1) or the ndim-dimensional IFFT (if isign = -1) of a complex vector data.
/// @param data Vector of data to be transformed. It has the length of nn[0] * nn[1] * nn[2] * ...  and will be replaced by the result at the end of the function. The complex data is stored as (real0, imag0, real1, imag1, ...). For a two-dimensional array, this is equivalent to storing the array by rows.
/// @param nn Vector of integers specifying the lengths of each dimension. All nn[i] must be powers of 2.
/// @param isign If isign is 1, compute the FFT, if isign is -1, compute the IFFT.
/// @param normalize The are three commonly used definitions of the Fourier transform:
//    If normalize = 0, the factor 1 / sqrt(n) is included in the forward and inverse transform.
//    If normalize = 1, the factor 1 / n is included in the forward transform (that is, the fourier coefficients get normalized).
//    If normalize = -1, the factor 1 / n is included in the inverse transform (that is, the inverse fourier coefficients get normalized).
// Note: for consistent results, the same value of normalize should be used for both the forward and inverse transform using the same data.
// The normalize factor used in numpy would be equivalent to set normalize = -1.
void fft_n(Vector<double> &data, const Vector<int> &nn, const int isign, const int normalize = 1);

// void fft_3_real(double *data, double *speq, const int isign, const int nn1, const int nn2, const int nn3);

// void fft_2_3_real(double *data, double *speq, const int isign, const int nn1, const int nn2, const int nn3);

// void fft_2_3_real(Matrix3d<double> &data, Matrix<double> &speq, const int isign);

// void fft_2_3_real(Matrix<double> &data, Vector<double> &speq, const int isign);

#endif  // FFT_HPP
