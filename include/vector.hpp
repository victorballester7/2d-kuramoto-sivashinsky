#ifndef VECTOR_HPP
#define VECTOR_HPP

#include <cmath>
#include <ostream>
#include <stdexcept>
#include <vector>

using namespace std;

// Note: const member functions cannot modify the object

void print() const;  // Print vector to stdout
int size() const;    // Return dimension of the vector
T norm() const;      // Return norm of the vector
template <class T>
vector<T> normalized() const;      // Return normalized vector
T dot(const vector<T>& v_) const;  // Return dot product with v_
template <class T>
vector<T> cross(const vector<T>& v_) const;  // Return cross product with v_

T operator[](int i) const;  // Access element i (const), i = 0, ..., n-1
T& operator[](int i);       // Access element i (non-const), i = 0, ..., n-1

template <class T>
vector<T> operator+(const vector<T>& v_) const;  // vector addition

template <class T>
vector<T> operator=(const vector<T> v_);  // Assignment

template <class T>
vector<T> operator+=(const vector<T>& v_);  // vector addition

template <class T>
vector<T> operator+(T a) const;  // Addition of a vector filled with a

template <class T>
vector<T> operator+=(T a) const;  // Addition of a vector filled with a

template <class T>
vector<T> operator-(const vector<T>& v_) const;  // vector subtraction

template <class T>
vector<T> operator-(T a) const;  // vector subtraction

template <class T>
vector<T> operator-=(const vector<T>& v_);  // vector subtraction

template <class T>
vector<T> operator*(T a) const;  // Scalar multiplication from the right

template <class T>
vector<T> operator*=(T a) const;  // Scalar multiplication from the right

template <class T>
// vector<T> operator*(T a, const vector<T>& v);    // Scalar multiplication from the left. We implement this as a non-member function
template <class T>
vector<T> operator*(const vector<T>& v_) const;  // vector multiplication (Hadamard product

template <class T>
vector<T> operator*=(const vector<T>& v_) const;  // vector multiplication (Hadamard product

template <class T>
vector<T> operator/(T a) const;  // Scalar division

template <class T>
vector<T> operator/=(T a) const;  // Scalar division

template <class T>
vector<T> operator/(const vector<T>& v_) const;  // vector division (Hadamard division)

template <class T>
vector<T> operator^(T a) const;  // Scalar power

template <class T>
bool operator==(const vector<T>& v_) const;  // vector equality

template <class T>
bool operator!=(const vector<T>& v_) const;  // vector inequality

template <class T>
friend std::ostream& operator<<(std::ostream& os, const vector<T>& Vec);  // Print vector to output stream

#endif  // LINEAR_ALGEBRA_HPP
