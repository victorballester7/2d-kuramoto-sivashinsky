#ifndef SCHEME_FD_HPP
#define SCHEME_FD_HPP

#include "misc.hpp"
#include "vector.hpp"

// coefficients for the IMEXRK4 scheme
const double alphaI[4] = {343038331393.0 / 1130875731271.0, 288176579239.0 / 1140253497719.0, 253330171251.0 / 677500478386.0, 189462239225.0 / 1091147436423.0};
const double betaI[4] = {35965327958.0 / 140127563663.0, 19632212512.0 / 2700543775099.0, -173747147147.0 / 351772688865.0, 91958533623.0 / 727726057489.0};
const double alphaE[4] = {14.0 / 25.0, 777974228744.0 / 1346157007247.0, 251277807242.0 / 1103637129625.0, 113091689455.0 / 220187950967.0};
const double betaE[4] = {0.0, -251352885992.0 / 790610919619.0, -383714262797.0 / 1103637129625.0, -403360439203.0 / 1888264787188.0};

// @brief: computes a step of the finite difference method used to solve the system
// @param x: the data to be updated. At the beginning it is the data at the previous time step and at the end it is the data at the next time step.
// @param t: time at which the data is written
// @param prm: parameters of the system
void stepFiniteDiff(vector<double> &x, double &t, const Args &prm);

void stepIMEXRK4(vector<double> &x, double &t, const Args &prm);
vector<double> N(vector<double> &x, const Args &prm);

#endif  // SCHEME_FD_HPP