#ifndef MISC_HPP
#define MISC_HPP

#include <fftw3.h>

#include <cmath>
#include <iostream>
#include <stdfloat>

#include "../include/energy.hpp"

using namespace std;

// @brief: computes the initial condition
// @param x: x-coordinate
// @param y: y-coordinate
// @return: the value of the initial condition at (x, y)
double u0(double x, double y);

double U(double x, double y, double t);
double g(double x, double y, double t);
// @brief: writes the data to a file
// @param x: the data to be written
// @param nn: number of points in each direction
// @param t: time at which the data is written
// @param file: file to write the data to
void write(double *x, int nx, int ny, double t, ofstream &file);

void write_bis(fftw_complex *x, int nx, int ny_complex, double t, ofstream &file);

// @brief: computes the initial condition on the whole vector
// @param x: the vector to be updated
// @param prm: parameters of the system
void set_data(double *x, int nx, int ny);

// @brief: computes the wave numbers
// @param kx: wave numbers in the x-direction: 0, 0, ..., 0, 1, 1, ..., 1, ..., nx/2, nx/2, ..., nx/2, -nx/2, -nx/2, ..., -nx/2, ..., -1, -1, ..., -1
// @param ky: wave numbers in the y-direction: 0, 1, ..., ny/2, 0, 1, ..., ny/2, ...
// @param nx: number of points in the x-direction
// @param ny_complex: number of points in the y-direction (we are not storing all the frequencies in the y-direction, because our transform is from real to complex and so we can make use of the symmetry in the coefficients and save memory) (ny_complex = ny / 2 + 1)
void set_wave_numbers(double *kx, double *ky, int nx, int ny_complex);

// @brief: computes the coefficients of the factor C multiplying each time the IMEX-Euler scheme
// @param C: coefficients of the factor C
// @param kx: wave numbers in the x-direction
// @param ky: wave numbers in the y-direction
// @param nx: number of points in the x-direction
// @param ny_complex: number of points in the y-direction (ny_complex = ny / 2 + 1)
// @param dt: time step
// @param nu1: parameter nu1 of the pde
// @param nu2: parameter nu2 of the pde
void set_C1(double *C1, double *kx, double *ky, int nx, int ny_complex, double dt, double nu1, double nu2, double c);

// @brief: computes the coefficients of the factor C multiplying each time the IMEX-Euler scheme
// @param C: coefficients of the factor C
// @param kx: wave numbers in the x-direction
// @param ky: wave numbers in the y-direction
// @param nx: number of points in the x-direction
// @param ny_complex: number of points in the y-direction (ny_complex = ny / 2 + 1)
// @param dt: time step
// @param nu1: parameter nu1 of the pde
// @param nu2: parameter nu2 of the pde
void set_C2(double *C2, double *kx, double *ky, int nx, int ny_complex, double dt, double nu1, double nu2, double c);

// @brief: computes the root of the Lagrange polynomial of degree 2 that passes through the points (x2 - 2 * dt, y0), (x2 - dt, y1), (x2, y2) assuming there is a root in the interval (x2 - 2 * dt, x2). Be aware that this function should not be used if you have not previously checked that there is a root in the interval (x2 - 2 * dt, x2).
// @param x2: x-coordinate of the last point
// @param dt: time step
// @param y0: value of the function at x2 - 2 * dt
// @param y1: value of the function at x2 - dt
// @param y2: value of the function at x2
// @return: the root of the Lagrange polynomial of degree 2 that passes through the points (x2 - 2 * dt, y0), (x2 - dt, y1), (x2, y2)
double lagrange_root(double x2, double dt, double y0, double y1, double y2);

// @brief: computes the value at x = t of the Lagrange polynomial of degree 2 that passes through the points (x2 - 2 * dt, y0), (x2 - dt, y1), (x2, y2).
// @param t: x-coordinate at which we want to evaluate the Lagrange polynomial
// @param x2: x-coordinate of the last point
// @param dt: time step
// @param y0: value of the function at x2 - 2 * dt
// @param y1: value of the function at x2 - dt
// @param y2: value of the function at x2
// @return: the value at x = t of the Lagrange polynomial of degree 2 that passes through the points (x2 - 2 * dt, y0), (x2 - dt, y1), (x2, y2)
double lagrange_eval(double t, double x2, double dt, double y0, double y1, double y2);

// @brief: computes the root of the Lagrange polynomial of degree 2 that passes through the points (x2 - 2 * dt, y0), (x2 - dt, y1), (x2, y2) assuming there is a root in the interval (x2 - 2 * dt, x2). Be aware that this function should not be used if you have not previously checked that there is a root in the interval (x2 - 2 * dt, x2).
// @param x2: x-coordinate of the last point
// @param dt: time step
// @param y0: value of the function at x2 - 2 * dt
// @param y1: value of the function at x2 - dt
// @param y2: value of the function at x2
// @return: the root of the Lagrange polynomial of degree 2 that passes through the points (x2 - 2 * dt, y0), (x2 - dt, y1), (x2, y2)
double lagrange_root(double x2, double dt, My_double y0, My_double y1, My_double y2);

// @brief: computes the value at x = t of the Lagrange polynomial of degree 2 that passes through the points (x2 - 2 * dt, y0), (x2 - dt, y1), (x2, y2).
// @param t: x-coordinate at which we want to evaluate the Lagrange polynomial
// @param x2: x-coordinate of the last point
// @param dt: time step
// @param y0: value of the function at x2 - 2 * dt
// @param y1: value of the function at x2 - dt
// @param y2: value of the function at x2
// @return: the value at x = t of the Lagrange polynomial of degree 2 that passes through the points (x2 - 2 * dt, y0), (x2 - dt, y1), (x2, y2)
double lagrange_eval(double t, double x2, double dt, My_double y0, My_double y1, My_double y2);
void print(string str, int64_t time);
#endif  // MISC_HPP