#include "../include/vector.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace std;

template <class T>
T norm(const vector<T>& v) {
  T norm = 0.;
  for (int i = 0; i < v.size(); i++)
    norm += v[i] * v[i];
  return sqrt(norm);
  // T squaredSum = transform_reduce(v.begin(), v.end(), v.begin(), T(0), std::plus<>(), [](const T& value) {
  //   return value * value;
  // });
  // return sqrt(squaredSum);
}

template <class T>
vector<T> normalized(const vector<T>& v) {
  T norm = norm(v);
  vector<T> u(v.size());
  if (norm == 0.)
    throw std::invalid_argument("normalized: norm is zero");
  for (int i = 0; i < v.size(); i++)
    u.v[i] = v[i] / norm;
  return u;
}

template <class T>
T operator|(const vector<T>& u, const vector<T>& v) {
  if (u.size() != v.size())
    throw std::invalid_argument("operator*: dimension mismatch");
  T dot = 0.;
  for (int i = 0; i < u.size(); i++)
    dot += u[i] * v[i];
  return dot;
}

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

template <class T>
vector<T> operator+(const vector<T>& u, const vector<T>& v) {
  if (u.size() != v.size())
    throw std::invalid_argument("operator+: dimension mismatch");
  vector<T> w(u.size());
  for (int i = 0; i < u.size(); i++)
    w[i] = u[i] + v[i];
  return w;
}

template <class T>
vector<T> operator+(const vector<T>& v, const T a) {
  vector<T> w(v.size());
  for (int i = 0; i < v.size(); i++) {
    w[i] = v[i] + a;
  }
  return w;
}

template <class T>
vector<T> operator+(const T a, const vector<T>& v) {
  return v + a;
}

template <class T>
vector<T> operator+=(vector<T>& u, const vector<T>& v) {
  if (u.size() != v.size())
    throw std::invalid_argument("operator+: dimension mismatch");
  for (int i = 0; i < u.size(); i++)
    u[i] += v[i];
  return u;
}

template <class T>
vector<T> operator+=(vector<T>& v, const T a) {
  for (int i = 0; i < v.size(); i++)
    v[i] += a;
  return v;
}

template <class T>
vector<T> operator-(const vector<T>& u, const vector<T>& v) {
  if (u.size() != v.size())
    throw std::invalid_argument("operator-: dimension mismatch");
  vector<T> w(u.size());
  for (int i = 0; i < u.size(); i++)
    w[i] = u[i] - v[i];
  return w;
}

template <class T>
vector<T> operator-(const vector<T>& v, const T a) {
  return v + (-a);
}

template <class T>
vector<T> operator-(const T a, const vector<T>& v) {
  return v + (-a);
}

template <class T>
vector<T> operator-=(vector<T>& u, const vector<T>& v) {
  if (u.size() != v.size())
    throw std::invalid_argument("operator-: dimension mismatch");
  for (int i = 0; i < u.size(); i++)
    u[i] -= v[i];
  return u;
}

template <class T>
vector<T> operator-=(vector<T>& v, const T a) {
  for (int i = 0; i < v.size(); i++)
    v[i] -= a;
  return v;
}

template <class T>
vector<T> operator*(const vector<T>& v, const T a) {
  vector<T> w(v.size());
  for (int i = 0; i < v.size(); i++)
    w[i] = v[i] * a;
  return w;
}

template <class T>
vector<T> operator*(const T a, const vector<T>& v) {
  return v * a;
}

template <class T>
vector<T> operator*=(vector<T>& v, const T a) {
  for (int i = 0; i < v.size(); i++)
    v[i] *= a;
  return v;
}

template <class T>
vector<T> operator*(const vector<T>& u, const vector<T>& v) {
  if (u.size() != v.size())
    throw std::invalid_argument("operator*: dimension mismatch");
  vector<T> w(u.size());
  for (int i = 0; i < u.size(); i++)
    w[i] = u[i] * v[i];
  return w;
}

template <class T>
vector<T> operator*=(vector<T>& u, const vector<T>& v) {
  if (u.size() != v.size())
    throw std::invalid_argument("operator*: dimension mismatch");
  for (int i = 0; i < u.size(); i++)
    u[i] *= v[i];
  return u;
}

template <class T>
vector<T> operator/(const vector<T>& v, const T a) {
  if (a == 0.)
    throw std::invalid_argument("operator/: division by zero");
  vector<T> w(v.size());
  for (int i = 0; i < v.size(); i++)
    w[i] = v[i] / a;
  return w;
}

template <class T>
vector<T> operator/(const T a, const vector<T>& v) {
  vector<T> w(v.size());
  for (int i = 0; i < v.size(); i++)
    w[i] = a / v[i];
  return w;
}

template <class T>
vector<T> operator/=(vector<T>& v, const T a) {
  if (a == 0.)
    throw std::invalid_argument("operator/: division by zero");
  for (int i = 0; i < v.size(); i++)
    v[i] /= a;
  return v;
}

template <class T>
vector<T> operator/(const vector<T>& u, const vector<T>& v) {
  if (u.size() != v.size())
    throw std::invalid_argument("operator*: dimension mismatch");
  vector<T> w(u.size());
  for (int i = 0; i < u.size(); i++) {
    if (v[i] == 0.)
      throw std::invalid_argument("operator/: division by zero");
    w[i] = u[i] / v[i];
  }
  return w;
}

template <class T>
vector<T> operator/=(vector<T>& u, const vector<T>& v) {
  if (u.size() != v.size())
    throw std::invalid_argument("operator*: dimension mismatch");
  for (int i = 0; i < u.size(); i++) {
    if (v[i] == 0.)
      throw std::invalid_argument("operator/: division by zero");
    u[i] /= v[i];
  }
  return u;
}

template <class T>
vector<T> operator^(const vector<T>& v, const T a) {
  vector<T> w(v.size());
  for (int i = 0; i < v.size(); i++)
    w[i] = pow(v[i], a);
  return w;
}

template <class T>
ostream& operator<<(std::ostream& os, const vector<T>& v) {
  os << "[";
  for (int i = 0; i < v.size() - 1; i++)
    os << v[i] << ", ";
  os << v[v.size() - 1] << "]";
  return os;
}