#ifndef LINEAR_ALGEBRA_HPP
#define LINEAR_ALGEBRA_HPP

#include <cmath>
#include <ostream>
#include <stdexcept>

template <class T>
class Vector {
 private:
  int n;  // Dimension
  T* v;   // Vector v_(n)

 public:
  // Constructors
  Vector() : n(0), v(NULL) {}  // Null vector (0-dim)
  Vector(int n_) : n(n_), v(new T[n_]{}) {}
  Vector(T* p, int n_) : n(n_), v(new T[n_]) {
    for (int i = 0; i < n_; i++) v[i] = p[i];
  }
  // template <typename... Args>
  // Vector(Args... args) : n(sizeof...(Args)), v(new T[n]{args...}) {}
  template <typename... Args>
  Vector(Args... args) : n(sizeof...(Args)), v(new T[n]{args...}) {}
  Vector(const Vector<T>& v_);  // Copy constructor

  // Destructors
  ~Vector();

  // Note: const member functions cannot modify the object
  void print() const;                          // Print vector to stdout
  int size() const;                            // Return dimension of the vector
  T norm() const;                              // Return norm of the vector
  Vector<T> normalized() const;                // Return normalized vector
  T dot(const Vector<T>& v_) const;            // Return dot product with v_
  Vector<T> cross(const Vector<T>& v_) const;  // Return cross product with v_

  T operator[](int i) const;                       // Access element i (const), i = 0, ..., n-1
  T& operator[](int i);                            // Access element i (non-const), i = 0, ..., n-1
  Vector<T> operator=(const Vector<T> v_);         // Assignment
  Vector<T> operator+(const Vector<T>& v_) const;  // Vector addition
  Vector<T> operator+=(const Vector<T>& v_);       // Vector addition
  Vector<T> operator+(T a) const;                  // Addition of a vector filled with a
  Vector<T> operator+=(T a) const;                 // Addition of a vector filled with a
  Vector<T> operator-(const Vector<T>& v_) const;  // Vector subtraction
  Vector<T> operator-(T a) const;                  // Vector subtraction
  Vector<T> operator-=(const Vector<T>& v_);       // Vector subtraction
  Vector<T> operator*(T a) const;                  // Scalar multiplication from the right
  Vector<T> operator*=(T a) const;                 // Scalar multiplication from the right
  // Vector<T> operator*(T a, const Vector<T>& v);    // Scalar multiplication from the left. We implement this as a non-member function
  Vector<T> operator*(const Vector<T>& v_) const;   // Vector multiplication (Hadamard product
  Vector<T> operator*=(const Vector<T>& v_) const;  // Vector multiplication (Hadamard product
  Vector<T> operator/(T a) const;                   // Scalar division
  Vector<T> operator/=(T a) const;                  // Scalar division
  Vector<T> operator/(const Vector<T>& v_) const;   // Vector division (Hadamard division)
  Vector<T> operator^(T a) const;                   // Scalar power
  bool operator==(const Vector<T>& v_) const;       // Vector equality
  bool operator!=(const Vector<T>& v_) const;       // Vector inequality
  template <class U>
  friend std::ostream& operator<<(std::ostream& os, const Vector<U>& Vec);  // Print vector to output stream
};

template <class T>
class Matrix {
 private:
  int m;  // First dimension (number of rows)
  int n;  // Second dimension (number of columns)
  T** M;  // Matrix M(m,n)

 public:
  // Constructors
  Matrix() : m(0), n(0), M(NULL) {}    // Null matrix (0 x 0)
  Matrix(int m_, int n_);              // Zero m x n matrix
  Matrix(const T* p, int m_, int n_);  // m x n matrix with elements p
  // Matrix(const T** p, int m_, int n_);  // m x n matrix with elements p
  Matrix(const Matrix<T>& M_);  // Copy constructor
  Matrix(int ax, T angle);      // rotation matrix around axis ax (ax = 1=x, 2=y, 3=z)

  // Destructor
  ~Matrix();

  // Assignment
  void print() const;           // Print matrix to stdout
  int nrows() const;            // Return number of rows
  int ncols() const;            // Return number of columns
  Matrix<T> identity(int n_);   // Return n x n identity matrix
  Matrix<T> transpose() const;  // Return transpose of the matrix
  // Matrix Inverse() const; // Return inverse of the matrix
  // double Determinant() const; // Return determinant of the matrix

  T operator()(int i, int j) const;  // Access element (i,j) (const)
  T& operator()(int i, int j);       // Access element (i,j) (non-const)

  Matrix<T>& operator=(const Matrix<T>& M_);  // Assignment

  Matrix<T> operator+(const Matrix<T>& M) const;  // Matrix addition
  Matrix<T> operator-(const Matrix<T>& M) const;  // Matrix subtraction
  Matrix<T> operator*(T a) const;                 // Scalar multiplication
  Matrix<T> operator/(T a) const;                 // Scalar division
  bool operator==(const Matrix<T>& M) const;      // Matrix equality
  bool operator!=(const Matrix<T>& M) const;      // Matrix inequality

  Matrix<T> operator*(const Matrix<T>& M_) const;  // Matrix multiplication
  Vector<T> operator*(const Vector<T>& v) const;   // Matrix-vector multiplication

  template <class U>
  friend std::ostream& operator<<(std::ostream& os, const Matrix<U>& Mat);  // Print matrix to output stream
};

template <class T>
Vector<T>::Vector(const Vector<T>& v_) : n(v_.n) {
  v = new T[n];
  for (int i = 0; i < n; i++) {
    v[i] = v_.v[i];
  }
}

template <class T>
Vector<T>::~Vector() {
  if (v != NULL) delete[] (v);
}

template <class T>
void Vector<T>::print() const {
  for (int i = 0; i < n; i++)
    printf("%f ", v[i]);
  printf("\n");
}

template <class T>
int Vector<T>::size() const {
  return n;
}

template <class T>
T Vector<T>::norm() const {
  T norm = 0.;
  for (int i = 0; i < n; i++)
    norm += v[i] * v[i];
  return sqrt(norm);
}

template <class T>
Vector<T> Vector<T>::normalized() const {
  T norm = this->norm();
  Vector<T> v_(n);
  for (int i = 0; i < n; i++)
    v_.v[i] = v[i] / norm;
  return v_;
}

template <class T>
T Vector<T>::dot(const Vector<T>& v_) const {
  if (n != v_.n)
    throw std::invalid_argument("Vector<T>::dot: dimension mismatch");
  T dot = 0.;
  for (int i = 0; i < n; i++)
    dot += v[i] * v_.v[i];
  return dot;
}

template <class T>
Vector<T> Vector<T>::cross(const Vector<T>& v_) const {
  if (n != 3 || v_.size() != 3)
    throw std::invalid_argument("Vector<T>::cross: dimension mismatch");
  Vector<T> u(3);
  u.v[0] = v[1] * v_(2) - v[2] * v_(1);
  u.v[1] = v[2] * v_(0) - v[0] * v_(2);
  u.v[2] = v[0] * v_(1) - v[1] * v_(0);
  return u;
}

template <class T>
T Vector<T>::operator[](int i) const {
  if (i < 0 || i >= n)
    throw std::out_of_range("Vector<T>::operator[]: index out of range");
  return v[i];
}

template <class T>
T& Vector<T>::operator[](int i) {
  if (i < 0 || i >= n)
    throw std::out_of_range("Vector<T>::operator[]: index out of range");
  return v[i];
}

template <class T>
Vector<T> Vector<T>::operator=(const Vector<T> v_) {
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
Vector<T> Vector<T>::operator+(const Vector<T>& v_) const {
  if (n != v_.n)
    throw std::invalid_argument("Vector<T>::operator+: dimension mismatch");
  Vector<T> u(n);
  for (int i = 0; i < n; i++)
    u.v[i] = v[i] + v_.v[i];
  return u;
}

template <class T>
Vector<T> Vector<T>::operator+=(const Vector<T>& v_) {
  if (n != v_.n)
    throw std::invalid_argument("Vector<T>::operator+: dimension mismatch");
  for (int i = 0; i < n; i++)
    v[i] += v_.v[i];
  return *this;
}

template <class T>
Vector<T> Vector<T>::operator+(T a) const {
  Vector<T> u(n);
  for (int i = 0; i < n; i++)
    u.v[i] = v[i] + a;
  return u;
}

template <class T>
Vector<T> operator+(T a, const Vector<T>& v) {
  Vector<T> result(v.size());
  for (int i = 0; i < v.size(); i++) {
    result[i] = a + v[i];
  }
  return result;
}

template <class T>
Vector<T> Vector<T>::operator+=(T a) const {
  for (int i = 0; i < n; i++)
    v[i] += a;
  return *this;
}

template <class T>
Vector<T> Vector<T>::operator-(const Vector<T>& v_) const {
  if (n != v_.n)
    throw std::invalid_argument("Vector<T>::operator-: dimension mismatch");
  Vector<T> u(n);
  for (int i = 0; i < n; i++)
    u.v[i] = v[i] - v_.v[i];
  return u;
}

template <class T>
Vector<T> Vector<T>::operator-(T a) const {
  Vector<T> u(n);
  for (int i = 0; i < n; i++)
    u.v[i] = v[i] - a;
  return u;
}

template <class T>
Vector<T> operator-(T a, const Vector<T>& v) {
  Vector<T> result(v.size());
  for (int i = 0; i < v.size(); i++) {
    result[i] = a - v[i];
  }
  return result;
}

template <class T>
Vector<T> Vector<T>::operator-=(const Vector<T>& v_) {
  if (n != v_.n)
    throw std::invalid_argument("Vector<T>::operator-: dimension mismatch");
  for (int i = 0; i < n; i++)
    v[i] -= v_.v[i];
  return *this;
}

template <class T>
Vector<T> Vector<T>::operator*(T a) const {
  Vector<T> u(n);
  for (int i = 0; i < n; i++)
    u.v[i] = v[i] * a;
  return u;
}

template <class T>
Vector<T> Vector<T>::operator*=(T a) const {
  for (int i = 0; i < n; i++)
    v[i] *= a;
  return *this;
}

template <class T>
Vector<T> operator*(T a, const Vector<T>& v) {
  Vector<T> result(v.size());
  for (int i = 0; i < v.size(); i++) {
    result[i] = a * v[i];
  }
  return result;
}

template <class T>
Vector<T> Vector<T>::operator*(const Vector<T>& v_) const {
  if (n != v_.n)
    throw std::invalid_argument("Vector<T>::operator*: dimension mismatch");
  Vector<T> u(n);
  for (int i = 0; i < n; i++)
    u[i] = v[i] * v_[i];
  return u;
}

template <class T>
Vector<T> Vector<T>::operator*=(const Vector<T>& v_) const {
  if (n != v_.n)
    throw std::invalid_argument("Vector<T>::operator*: dimension mismatch");
  for (int i = 0; i < n; i++)
    v[i] *= v_.v[i];
  return *this;
}

template <class T>
Vector<T> Vector<T>::operator/(T a) const {
  Vector<T> u(n);
  for (int i = 0; i < n; i++)
    u.v[i] = v[i] / a;
  return u;
}

template <class T>
Vector<T> operator/(T a, const Vector<T>& v) {
  Vector<T> result(v.size());
  for (int i = 0; i < v.size(); i++) {
    result[i] = a / v[i];
  }
  return result;
}

template <class T>
Vector<T> Vector<T>::operator/=(T a) const {
  for (int i = 0; i < n; i++)
    v[i] /= a;
  return *this;
}

template <class T>
Vector<T> Vector<T>::operator/(const Vector<T>& v_) const {
  if (n != v_.n)
    throw std::invalid_argument("Vector<T>::operator/: dimension mismatch");
  Vector<T> u(n);
  for (int i = 0; i < n; i++)
    u[i] = v[i] / v_[i];
  return u;
}

template <class T>
Vector<T> Vector<T>::operator^(T a) const {
  Vector<T> u(n);
  for (int i = 0; i < n; i++)
    u[i] = pow(v[i], a);
  return u;
}

template <class T>
bool Vector<T>::operator==(const Vector<T>& v_) const {
  if (n != v_.n)
    return false;
  for (int i = 0; i < n; i++) {
    if (v[i] != v_.v[i])
      return false;
  }
  return true;
}

template <class T>
bool Vector<T>::operator!=(const Vector<T>& v_) const {
  return !(*this == v_);
}

template <class T>
std::ostream& operator<<(std::ostream& os, const Vector<T>& u) {
  os << "(" << u(0) << ", " << u(1) << ", " << u(2) << ")";
  return os;
}

template <class T>
Matrix<T>::Matrix(int m_, int n_) : m(m_), n(n_) {
  // Allocate memory
  M = new T*[m_];
  for (int i = 0; i < m_; i++)
    M[i] = new T[n_];

  // Initialization
  for (int i = 0; i < m_; i++) {
    for (int j = 0; j < n_; j++)
      M[i][j] = 0.;
  }
}

template <class T>
Matrix<T>::Matrix(const T* p, int m_, int n_) : m(m_), n(n_) {
  // Allocate memory
  M = new T*[m_];
  for (int i = 0; i < m_; i++)
    M[i] = new T[n_];

  // Initialization
  for (int i = 0; i < m_; i++) {
    for (int j = 0; j < n_; j++)
      M[i][j] = p[i * n_ + j];
  }
}

template <class T>
Matrix<T>::Matrix(const Matrix<T>& M_) : m(M_.m), n(M_.n) {
  // Allocate memory
  M = new T*[m];
  for (int i = 0; i < m; i++)
    M[i] = new T[n];

  // Initialization
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++)
      M[i][j] = M_.M[i][j];
  }
}

template <class T>
Matrix<T>::Matrix(int ax, T angle) : m(3), n(3) {
  // Allocate memory
  M = new T*[m];
  for (int i = 0; i < m; i++)
    M[i] = new T[n];

  // Initialization
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++)
      M[i][j] = 0.;
  }

  // Rotation matrix around axis i
  if (ax == 1) {
    M[0][0] = 1.;
    M[1][1] = cos(angle);
    M[1][2] = sin(angle);
    M[2][1] = -sin(angle);
    M[2][2] = cos(angle);
  } else if (ax == 2) {
    M[0][0] = cos(angle);
    M[0][2] = -sin(angle);
    M[1][1] = 1.;
    M[2][0] = sin(angle);
    M[2][2] = cos(angle);
  } else if (ax == 3) {
    M[0][0] = cos(angle);
    M[0][1] = sin(angle);
    M[1][0] = -sin(angle);
    M[1][1] = cos(angle);
    M[2][2] = 1.;
  } else {
    throw std::invalid_argument("Invalid axis index");
  }
}

template <class T>
Matrix<T>::~Matrix() {
  for (int i = 0; i < m; i++)
    delete[] M[i];
  delete[] M;
}

template <class T>
void Matrix<T>::print() const {
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++)
      printf("%f ", M[i][j]);
    printf("\n");
  }
}

template <class T>
int Matrix<T>::nrows() const {
  return m;
}

template <class T>
int Matrix<T>::ncols() const {
  return n;
}

template <class T>
Matrix<T> Matrix<T>::identity(int n_) {
  Matrix<T> I(n_, n_);
  for (int i = 0; i < n_; i++)
    I.M[i][i] = 1.;
  return I;
}

template <class T>
Matrix<T> Matrix<T>::transpose() const {
  Matrix<T> result(n, m);
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++)
      result.M[j][i] = M[i][j];
  }
  return result;
}

template <class T>
T Matrix<T>::operator()(int i, int j) const {
  if (i < 0 || i >= m || j < 0 || j >= n)
    throw std::out_of_range("Matrix<T>::operator(): index out of range");
  return M[i][j];
}

template <class T>
T& Matrix<T>::operator()(int i, int j) {
  if (i < 0 || i >= m || j < 0 || j >= n)
    throw std::out_of_range("Matrix<T>::operator(): index out of range");
  return M[i][j];
}

template <class T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T>& M_) {
  if (m != M_.m || n != M_.n) {
    for (int i = 0; i < m; i++)
      delete[] M[i];
    delete[] M;
    m = M_.m;
    n = M_.n;
    M = new T*[m];
    for (int i = 0; i < m; i++)
      M[i] = new T[n];
  }
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++)
      M[i][j] = M_.M[i][j];
  }
  return *this;
}

template <class T>
Matrix<T> Matrix<T>::operator+(const Matrix<T>& M_) const {
  if (m != M_.m || n != M_.n)
    throw std::invalid_argument("Matrix<T>::operator+: incompatible dimensions");

  Matrix<T> result(m, n);
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++)
      result.M[i][j] = M[i][j] + M_.M[i][j];
  }
  return result;
}

template <class T>
Matrix<T> Matrix<T>::operator-(const Matrix<T>& M_) const {
  if (m != M_.m || n != M_.n)
    throw std::invalid_argument("Matrix<T>::operator-: incompatible dimensions");

  Matrix<T> result(m, n);
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++)
      result.M[i][j] = M[i][j] - M_.M[i][j];
  }
  return result;
}

template <class T>
Matrix<T> Matrix<T>::operator*(T a) const {
  Matrix<T> result(m, n);
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++)
      result.M[i][j] = M[i][j] * a;
  }
  return result;
}

template <class T>
Matrix<T> Matrix<T>::operator/(T a) const {
  Matrix<T> result(m, n);
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++)
      result.M[i][j] = M[i][j] / a;
  }
  return result;
}

template <class T>
Matrix<T> Matrix<T>::operator*(const Matrix<T>& M_) const {
  if (n != M_.m)
    throw std::invalid_argument("Matrix<T>::operator*: incompatible dimensions");

  Matrix<T> result(m, M_.n);
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < M_.n; j++) {
      for (int k = 0; k < n; k++)
        result.M[i][j] += M[i][k] * M_.M[k][j];
    }
  }
  return result;
}

template <class T>
Vector<T> Matrix<T>::operator*(const Vector<T>& v_) const {
  if (n != v_.size())
    throw std::invalid_argument("Matrix<T>::operator*: incompatible dimensions");

  Vector<T> result(m);
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++)
      result[i] += M[i][j] * v_(j);
  }
  return result;
}

template <class T>
std::ostream& operator<<(std::ostream& os, const Matrix<T>& M) {
  for (int i = 0; i < M.m; i++) {
    for (int j = 0; j < M.n; j++)
      os << M.M[i][j] << " ";
    os << std::endl;
  }
  return os;
}

#endif  // LINEAR_ALGEBRA_HPP
