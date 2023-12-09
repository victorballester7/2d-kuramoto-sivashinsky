#ifndef VECTOR_HPP
#define VECTOR_HPP

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace std;

// ------------------ NORM FUNCTIONS ------------------

// @brief Computes the norm of a vector.
// @param v Vector.
// @return Norm of the vector.
template <class T>
T norm(const vector<T>& v) {
  T norm = 0.;
  for (uint i = 0; i < v.size(); i++)
    norm += v[i] * v[i];
  return sqrt(norm);
  // T squaredSum = transform_reduce(v.begin(), v.end(), v.begin(), T(0), std::plus<>(), [](const T& value) {
  //   return value * value;
  // });
  // return sqrt(squaredSum);
}

// @brief Computes the normalized vector.
// @param v Vector.
// @return Normalized vector.
template <class T>
vector<T> normalized(const vector<T>& v) {
  T norm_v = norm(v);
  vector<T> u(v.size());
  if (norm_v == 0.)
    throw std::invalid_argument("normalized: norm is zero");
  for (uint i = 0; i < v.size(); i++)
    u[i] = v[i] / norm_v;
  return u;
}

// ----------------------------------------------------

// ------------------ MISCELLANEOUS -------------------

// @brief Computes the dot product of two vectors.
// @param u Vector.
// @param v Vector.
// @return Dot product of the two vectors.
template <class T>
T operator|(const vector<T>& u, const vector<T>& v) {
  if (u.size() != v.size())
    throw std::invalid_argument("operator*: dimension mismatch");
  T dot = 0.;
  for (uint i = 0; i < u.size(); i++)
    dot += u[i] * v[i];
  return dot;
}

// @brief Computes the cross product of two vectors.
// @param u Vector.
// @param v Vector.
// @return Cross product of the two vectors.
template <class T>
vector<T> operator^(const vector<T>& u, const vector<T>& v) {
  if (u.size() != 3 || v.size() != 3)
    throw std::invalid_argument("operator^: dimension mismatch");
  vector<T> w(3);
  w[0] = u[1] * v[2] - u[2] * v[1];
  w[1] = u[2] * v[0] - u[0] * v[2];
  w[2] = u[0] * v[1] - u[1] * v[0];
  return w;
}

// ----------------------------------------------------

// ------------------ ADDITION ------------------------

// @brief Computes the sum of two vectors.
// @param u Vector.
// @param v Vector.
// @return Sum of the two vectors.
template <class T>
vector<T> operator+(const vector<T>& u, const vector<T>& v) {
  if (u.size() != v.size())
    throw std::invalid_argument("operator+: dimension mismatch");
  vector<T> w(u.size());
  for (uint i = 0; i < u.size(); i++)
    w[i] = u[i] + v[i];
  return w;
}

// @brief Computes the sum of a vector and a scalar (Matlab-like syntax).
// @param v Vector.
// @param a Scalar.
// @return Sum of the vector and the scalar (vector of the same dimension with all elements equal to the scalar)
template <class T>
vector<T> operator+(const vector<T>& v, const T a) {
  vector<T> w(v.size());
  for (uint i = 0; i < v.size(); i++) {
    w[i] = v[i] + a;
  }
  return w;
}

// @brief Computes the sum of a scalar and a vector (Matlab-like syntax).
// @param a Scalar.
// @param v Vector.
// @return Sum of the scalar and the vector (vector of the same dimension with all elements equal to the scalar)
template <class T>
vector<T> operator+(const T a, const vector<T>& v) {
  return v + a;
}

// @brief Computes the sum of two vectors and stores the result in the first one.
// @param u Vector.
// @param v Vector.
// @return Sum of the two vectors.
template <class T>
vector<T> operator+=(vector<T>& u, const vector<T>& v) {
  if (u.size() != v.size())
    throw std::invalid_argument("operator+: dimension mismatch");
  for (uint i = 0; i < u.size(); i++)
    u[i] += v[i];
  return u;
}

// @brief Computes the sum of a vector and a scalar (Matlab-like syntax) and stores the result in the first one.
// @param v Vector.
// @param a Scalar.
// @return Sum of the vector and the scalar (vector of the same dimension with all elements equal to the scalar)
template <class T>
vector<T> operator+=(vector<T>& v, const T a) {
  for (uint i = 0; i < v.size(); i++)
    v[i] += a;
  return v;
}

// ----------------------------------------------------

// ------------------ SUBTRACTION ---------------------

// @brief Computes the difference of two vectors.
// @param u Vector.
// @param v Vector.
// @return Difference of the two vectors.
template <class T>
vector<T> operator-(const vector<T>& u, const vector<T>& v) {
  if (u.size() != v.size())
    throw std::invalid_argument("operator-: dimension mismatch");
  vector<T> w(u.size());
  for (uint i = 0; i < u.size(); i++)
    w[i] = u[i] - v[i];
  return w;
}

// @brief Computes the difference of a vector and a scalar (Matlab-like syntax).
// @param v Vector.
// @param a Scalar.
// @return Difference of the vector and the scalar (vector of the same dimension with all elements equal to the scalar)
template <class T>
vector<T> operator-(const vector<T>& v, const T a) {
  vector<T> w(v.size());
  for (uint i = 0; i < v.size(); i++)
    w[i] = v[i] - a;
  return w;
}

// @brief Computes the difference of a scalar and a vector (Matlab-like syntax).
// @param a Scalar.
// @param v Vector.
// @return Difference of the scalar and the vector (vector of the same dimension with all elements equal to the scalar)
template <class T>
vector<T> operator-(const T a, const vector<T>& v) {
  vector<T> w(v.size());
  for (uint i = 0; i < v.size(); i++)
    w[i] = a - v[i];
  return w;
}

// @brief Computes the difference of two vectors and stores the result in the first one.
// @param u Vector.
// @param v Vector.
// @return Difference of the two vectors.
template <class T>
vector<T> operator-=(vector<T>& u, const vector<T>& v) {
  if (u.size() != v.size())
    throw std::invalid_argument("operator-: dimension mismatch");
  for (uint i = 0; i < u.size(); i++)
    u[i] -= v[i];
  return u;
}

// @brief Computes the difference of a vector and a scalar (Matlab-like syntax) and stores the result in the first one.
// @param v Vector.
// @param a Scalar.
// @return Difference of the vector and the scalar (vector of the same dimension with all elements equal to the scalar)
template <class T>
vector<T> operator-=(vector<T>& v, const T a) {
  for (uint i = 0; i < v.size(); i++)
    v[i] -= a;
  return v;
}

// ----------------------------------------------------

// ------------------ MULTIPLICATION ------------------

// @brief Computes the product of two vectors (element-wise, Matlab-like syntax).
// @param u Vector.
// @param v Vector.
// @return Product of the two vectors.
template <class T>
vector<T> operator*(const vector<T>& v, const T a) {
  vector<T> w(v.size());
  for (uint i = 0; i < v.size(); i++)
    w[i] = v[i] * a;
  return w;
}

// @brief Computes the product of a vector and a scalar (Matlab-like syntax).
// @param v Vector.
// @param a Scalar.
// @return Product of the vector and the scalar (vector of the same dimension with all elements equal to the scalar)
template <class T>
vector<T> operator*(const T a, const vector<T>& v) {
  return v * a;
}

// @brief Computes the product of two vectors (element-wise, Matlab-like syntax) and stores the result in the first one.
// @param u Vector.
// @param v Vector.
// @return Product of the two vectors.
template <class T>
vector<T> operator*=(vector<T>& v, const T a) {
  for (uint i = 0; i < v.size(); i++)
    v[i] *= a;
  return v;
}

// @brief Computes the product of two vectors (element-wise, Matlab-like syntax).
// @param u Vector.
// @param v Vector.
// @return Product of the two vectors.
template <class T>
vector<T> operator*(const vector<T>& u, const vector<T>& v) {
  if (u.size() != v.size())
    throw std::invalid_argument("operator*: dimension mismatch");
  vector<T> w(u.size());
  for (uint i = 0; i < u.size(); i++)
    w[i] = u[i] * v[i];
  return w;
}

// @brief Computes the product of two vectors (element-wise, Matlab-like syntax) and stores the result in the first one.
// @param u Vector.
// @param v Vector.
// @return Product of the two vectors.
template <class T>
vector<T> operator*=(vector<T>& u, const vector<T>& v) {
  if (u.size() != v.size())
    throw std::invalid_argument("operator*: dimension mismatch");
  for (uint i = 0; i < u.size(); i++)
    u[i] *= v[i];
  return u;
}

// ----------------------------------------------------

// ------------------ DIVISION ------------------------

// @brief Computes the division of two vectors (element-wise, Matlab-like syntax).
// @param u Vector.
// @param v Vector.
// @return Division of the two vectors.
template <class T>
vector<T> operator/(const vector<T>& v, const T a) {
  if (a == 0.)
    throw std::invalid_argument("operator/: division by zero");
  vector<T> w(v.size());
  for (uint i = 0; i < v.size(); i++)
    w[i] = v[i] / a;
  return w;
}

// @brief Computes the division of a vector and a scalar (Matlab-like syntax).
// @param v Vector.
// @param a Scalar.
// @return Division of the vector and the scalar (vector of the same dimension with all elements equal to the scalar)
template <class T>
vector<T> operator/(const T a, const vector<T>& v) {
  vector<T> w(v.size());
  for (uint i = 0; i < v.size(); i++)
    w[i] = a / v[i];
  return w;
}

// @brief Computes the division of two vectors (element-wise, Matlab-like syntax) and stores the result in the first one.
// @param u Vector.
// @param v Vector.
// @return Division of the two vectors.
template <class T>
vector<T> operator/=(vector<T>& v, const T a) {
  if (a == 0.)
    throw std::invalid_argument("operator/: division by zero");
  for (uint i = 0; i < v.size(); i++)
    v[i] /= a;
  return v;
}

// @brief Computes the division of two vectors (element-wise, Matlab-like syntax).
// @param u Vector.
// @param v Vector.
// @return Division of the two vectors.
template <class T>
vector<T> operator/(const vector<T>& u, const vector<T>& v) {
  if (u.size() != v.size())
    throw std::invalid_argument("operator*: dimension mismatch");
  vector<T> w(u.size());
  for (uint i = 0; i < u.size(); i++) {
    if (v[i] == 0.)
      throw std::invalid_argument("operator/: division by zero");
    w[i] = u[i] / v[i];
  }
  return w;
}

// @brief Computes the division of two vectors (element-wise, Matlab-like syntax) and stores the result in the first one.
// @param u Vector.
// @param v Vector.
// @return Division of the two vectors.
template <class T>
vector<T> operator/=(vector<T>& u, const vector<T>& v) {
  if (u.size() != v.size())
    throw std::invalid_argument("operator*: dimension mismatch");
  for (uint i = 0; i < u.size(); i++) {
    if (v[i] == 0.)
      throw std::invalid_argument("operator/: division by zero");
    u[i] /= v[i];
  }
  return u;
}

// ----------------------------------------------------

// ------------------ POWER ---------------------------

// @brief Computes the power of a vector (element-wise, Matlab-like syntax).
// @param v Vector.
// @param a Scalar.
// @return Power of the vector (vector of the same dimension with all elements equal to the scalar)
template <class T>
vector<T> operator^(const vector<T>& v, const T a) {
  vector<T> w(v.size());
  for (uint i = 0; i < v.size(); i++)
    w[i] = pow(v[i], a);
  return w;
}

// ----------------------------------------------------

// ------------------ PRINT ---------------------------

// @brief Prints a vector.
// @param os Output stream.
// @param v Vector.
// @return Output stream.
template <class T>
ostream& operator<<(std::ostream& os, const vector<T>& v) {
  os << "[";
  for (uint i = 0; i < v.size() - 1; i++)
    os << v[i] << ", ";
  os << v[v.size() - 1] << "]";
  return os;
}

#endif  // VECTOR_HPP
