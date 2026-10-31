// Name: Rémi, Roux
// Date: 31 October 2026 (35mins)
//
// Code output:
// **************************************************************
// xlo = 0.0000000e+00; xhi = 3.5000000e+00; nint = 7; nfv = 21
// xlo = 0.0000000e+00; xhi = 3.2640832e+00; nint = 7; nfv = 42
// xlo = 0.0000000e+00; xhi = 3.2634687e+00; nint = 7; nfv = 63
// xlo = 0.0000000e+00; xhi = 3.2634687e+00; nint = 7; nfv = 84
// Particle position (t = 5.0000000e+00): 3.2634687e+00
// **************************************************************

#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

#include "../../libs/integrators.h"
#include "../../libs/root_finders.h"

using namespace std;

double g_T = 5;

double Time(double);
double TimeMinusGT(double);
double Inv_Velocity(double);

int main() {
  double root;
  // by studying TimeMinusGT against x we see that 
  // the root is between x=3 and x=4
  double a = 3;
  double b = 4;

  double xtol = 1e-6;
  double ytol = -1;

  cout << setiosflags(ios_base::scientific) << setprecision(7);

  //fdata.open("trajectory_inversion.dat");
  
  //double y = 0;
  //for (int i = 0; i < 1000; i++) {
  //  y = i * 0.01;
  //  fdata << y << "  " << TimeMinusGT(y) << endl;
  //}

  //fdata.close();

  // Inv_Velocity is the derivative of TimeMinus (because of its definition
  // through an integral)
  if (Newton(TimeMinusGT, Inv_Velocity, a, b, xtol, ytol, root)) {
    cout << "some error occurred in Newton integration" << endl;
  }

  cout << "Particle position (t = " << g_T << "): " << root << endl;

  return 0;
}

double Time(double x) {
  // extremes of integration in the defintion of t(x)
  double xlo = 0.0;
  // double xhi = x    (why waste any memory?)

  int Ng   = 3;
  int nint = ceil(fabs(x - xlo) * 2.0);

  static int nfv = 0; // cumulative number of function evalutations   

  double time = GaussQuad(Inv_Velocity, xlo, x, nint, Ng);     
  nfv += Ng * nint;   // Inv_Velocity gets called N_g * nint times
  
  cout << "xlo = " << xlo << "; xhi = " << x 
       << "; nint = " << nint  << "; nfv = " << nfv << endl;   

  return time;
}

// computes Time - g_T: we want to find particles position at T = g_T
// thus we want to find the root of the function T - g_T
double TimeMinusGT(double x) {
  return Time(x) - g_T;
}

// computes inverse of particle's velocity
double Inv_Velocity(double x) {
  return sqrt(sin(x) * sin(x) * sin(x) * sin(x) + 2.0);
}
