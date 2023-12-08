#include "../include/vector.hpp"

#include <iostream>

template <class T>
vector<T>::vector(const vector<T>& v_) : n(v_.n) {
  v = new T[n];
  for (int i = 0; i < n; i++) {
    v[i] = v_.v[i];
  }
}

template <class T>
vector<T>::~vector() {
  if (v != NULL) delete[] (v);
}

template <class T>
void vector<T>::print() const {
  for (int i = 0; i < n; i++)
    printf("%f ", v[i]);
  printf("\n");
}

template <class T>
int vector<T>::size() const {
  return n;
}

template <class T>
T vector<T>::norm() const {
  T norm = 0.;
  for (int i = 0; i < n; i++)
    norm += v[i] * v[i];
  return sqrt(norm);
}

template <class T>
vector<T> vector<T>::normalized() const {
  T norm = this->norm();
  vector<T> v_(n);
  for (int i = 0; i < n; i++)
    v_.v[i] = v[i] / norm;
  return v_;
}

template <class T>
T vector<T>::dot(const vector<T>& v_) const {
  if (n != v_.n)
    throw std::invalid_argument("vector<T>::dot: dimension mismatch");
  T dot = 0.;
  for (int i = 0; i < n; i++)
    dot += v[i] * v_.v[i];
  return dot;
}

template <class T>
vector<T> vector<T>::cross(const vector<T>& v_) const {
  if (n != 3 || v_.size() != 3)
    throw std::invalid_argument("vector<T>::cross: dimension mismatch");
  vector<T> u(3);
  u.v[0] = v[1] * v_(2) - v[2] * v_(1);
  u.v[1] = v[2] * v_(0) - v[0] * v_(2);
  u.v[2] = v[0] * v_(1) - v[1] * v_(0);
  return u;
}

template <class T>
T vector<T>::operator[](int i) const {
  if (i < 0 || i >= n)
    throw std::out_of_range("vector<T>::operator[]: index out of range");
  return v[i];
}

template <class T>
T& vector<T>::operator[](int i) {
  if (i < 0 || i >= n)
    throw std::out_of_range("vector<T>::operator[]: index out of range");
  return v[i];
}

template <class T>
vector<T> vector<T>::operator=(const vector<T> v_) {
  if (n != v_.n) {
    if (v != NULL) delete[] (v);
    n = v_.n;
    v = new T[n];
  }
  for (int i = 0; i < n; i++)
    v[i] = v_.v[i];
  return *this;
}

template <class T>
vector<T> vector<T>::operator+(const vector<T>& v_) const {
  if (n != v_.n)
    throw std::invalid_argument("vector<T>::operator+: dimension mismatch");
  vector<T> u(n);
  for (int i = 0; i < n; i++)
    u.v[i] = v[i] + v_.v[i];
  return u;
}

template <class T>
vector<T> vector<T>::operator+=(const vector<T>& v_) {
  if (n != v_.n)
    throw std::invalid_argument("vector<T>::operator+: dimension mismatch");
  for (int i = 0; i < n; i++)
    v[i] += v_.v[i];
  return *this;
}

template <class T>
vector<T> vector<T>::operator+(T a) const {
  vector<T> u(n);
  for (int i = 0; i < n; i++)
    u.v[i] = v[i] + a;
  return u;
}

template <class T>
vector<T> operator+(T a, const vector<T>& v) {
  vector<T> result(v.size());
  for (int i = 0; i < v.size(); i++) {
    result[i] = a + v[i];
  }
  return result;
}

template <class T>
vector<T> vector<T>::operator+=(T a) const {
  for (int i = 0; i < n; i++)
    v[i] += a;
  return *this;
}

template <class T>
vector<T> vector<T>::operator-(const vector<T>& v_) const {
  if (n != v_.n)
    throw std::invalid_argument("vector<T>::operator-: dimension mismatch");
  vector<T> u(n);
  for (int i = 0; i < n; i++)
    u.v[i] = v[i] - v_.v[i];
  return u;
}

template <class T>
vector<T> vector<T>::operator-(T a) const {
  vector<T> u(n);
  for (int i = 0; i < n; i++)
    u.v[i] = v[i] - a;
  return u;
}

template <class T>
vector<T> operator-(T a, const vector<T>& v) {
  vector<T> result(v.size());
  for (int i = 0; i < v.size(); i++) {
    result[i] = a - v[i];
  }
  return result;
}

template <class T>
vector<T> vector<T>::operator-=(const vector<T>& v_) {
  if (n != v_.n)
    throw std::invalid_argument("vector<T>::operator-: dimension mismatch");
  for (int i = 0; i < n; i++)
    v[i] -= v_.v[i];
  return *this;
}

template <class T>
vector<T> vector<T>::operator*(T a) const {
  vector<T> u(n);
  for (int i = 0; i < n; i++)
    u.v[i] = v[i] * a;
  return u;
}

template <class T>
vector<T> vector<T>::operator*=(T a) const {
  for (int i = 0; i < n; i++)
    v[i] *= a;
  return *this;
}

template <class T>
vector<T> operator*(T a, const vector<T>& v) {
  vector<T> result(v.size());
  for (int i = 0; i < v.size(); i++) {
    result[i] = a * v[i];
  }
  return result;
}

template <class T>
vector<T> vector<T>::operator*(const vector<T>& v_) const {
  if (n != v_.n)
    throw std::invalid_argument("vector<T>::operator*: dimension mismatch");
  vector<T> u(n);
  for (int i = 0; i < n; i++)
    u[i] = v[i] * v_[i];
  return u;
}

template <class T>
vector<T> vector<T>::operator*=(const vector<T>& v_) const {
  if (n != v_.n)
    throw std::invalid_argument("vector<T>::operator*: dimension mismatch");
  for (int i = 0; i < n; i++)
    v[i] *= v_.v[i];
  return *this;
}

template <class T>
vector<T> vector<T>::operator/(T a) const {
  vector<T> u(n);
  for (int i = 0; i < n; i++)
    u.v[i] = v[i] / a;
  return u;
}

template <class T>
vector<T> operator/(T a, const vector<T>& v) {
  vector<T> result(v.size());
  for (int i = 0; i < v.size(); i++) {
    result[i] = a / v[i];
  }
  return result;
}

template <class T>
vector<T> vector<T>::operator/=(T a) const {
  for (int i = 0; i < n; i++)
    v[i] /= a;
  return *this;
}

template <class T>
vector<T> vector<T>::operator/(const vector<T>& v_) const {
  if (n != v_.n)
    throw std::invalid_argument("vector<T>::operator/: dimension mismatch");
  vector<T> u(n);
  for (int i = 0; i < n; i++)
    u[i] = v[i] / v_[i];
  return u;
}

template <class T>
vector<T> vector<T>::operator^(T a) const {
  vector<T> u(n);
  for (int i = 0; i < n; i++)
    u[i] = pow(v[i], a);
  return u;
}

template <class T>
bool vector<T>::operator==(const vector<T>& v_) const {
  if (n != v_.n)
    return false;
  for (int i = 0; i < n; i++) {
    if (v[i] != v_.v[i])
      return false;
  }
  return true;
}

template <class T>
bool vector<T>::operator!=(const vector<T>& v_) const {
  return !(*this == v_);
}

template <class T>
std::ostream& operator<<(std::ostream& os, const vector<T>& u) {
  os << "(" << u(0) << ", " << u(1) << ", " << u(2) << ")";
  return os;
}