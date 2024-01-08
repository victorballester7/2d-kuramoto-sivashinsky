#ifndef MISC_HPP
#define MISC_HPP

#include <fftw3.h>

#include <cmath>
#include <iostream>

// #include "fftw3.h"
using namespace std;

// @brief: computes the initial condition
// @param x: x-coordinate
// @param y: y-coordinate
// @return: the value of the initial condition at (x, y)
double u0(double x, double y);

// @brief: writes the data to a file
// @param x: the data to be written
// @param nn: number of points in each direction
// @param t: time at which the data is written
// @param file: file to write the data to
// @param print_complex: whether to print the imaginary part of the data or not
void write(double *x, int nx, int ny, double t, ofstream &file);

// void write_old(fftw_complex *x, int nx, int ny, double t, ofstream &file);
// void write_old_old(vector<double> &x, vector<uint> nn, double t, ofstream &file, bool print_complex);
// @brief: computes the mean of the data as the double integral of the data over the domain divided by the area of the domain
// @param x: the data
// @param prm: parameters of the system
// @return: the mean
// double mean(const vector<double> &x, const Args &prm);
double mean(double *x, int nx, int ny);
// @brief: computes the initial condition on the whole vector
// @param x: the vector to be updated
// @param prm: parameters of the system
void set_data(double *x, int nx, int ny);

void set_wave_numbers(double *kx, double *ky, int nx, int ny_complex);

void set_C(double *tmp, double *kx, double *ky, int nx, int ny_complex, double dt, double nu1, double nu2);

double root_dE_Lagrange(double x2, double dt, double y0, double y1, double y2);

double E_Lagrange(double t, double x2, double dt, double y0, double y1, double y2);
#endif  // MISC_HPP