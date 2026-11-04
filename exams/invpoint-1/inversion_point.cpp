// Name: Rémi Roux
// Date: 4 November 2026 (20 mins)
//
// Code Output:
// *****************************************************
// tlo = 0.0000000e+00; thi = 2.5000000e+00; nint = 5; Func eval = 20
// tlo = 2.5000000e+00; thi = 4.3197123e+00; nint = 4; Func eval = 36
// tlo = 4.3197123e+00; thi = 6.2836519e+00; nint = 4; Func eval = 52
// tlo = 6.2836519e+00; thi = 7.5449171e+00; nint = 3; Func eval = 64
// tlo = 7.5449171e+00; thi = 7.8471410e+00; nint = 1; Func eval = 68
// tlo = 7.8471410e+00; thi = 7.8597039e+00; nint = 1; Func eval = 72
// tlo = 7.8597039e+00; thi = 7.8597240e+00; nint = 1; Func eval = 76
// tlo = 7.8597240e+00; thi = 7.8597240e+00; nint = 1; Func eval = 80
// Inversion time (zero) = 7.8597240e+00
// *****************************************************

#include <iostream>
#include <iomanip>
#include <cmath>

#include "../../libs/integrators.h"
#include "../../libs/root_finders.h"

using namespace std;

double g_starttime = 0;
double g_startvel  = -1;

double Velocity(double);
double Acceleration(double);

int main() {
  double root;
  double a = 0;
  double b = 5;
  double xtol = 1e-6;
  double ytol = -1;

  // setting number output format
  cout << setiosflags(ios_base::scientific) << setprecision(7);

  Newton(Velocity, Acceleration, a, b, xtol, ytol, root);

  cout << "Inversion time (zero) = " << root << endl;

  return 0;
}

double Velocity(double t){ 
  static int nfv = 0; // Cumulative number of function evaluations ...  

  double tlo = g_starttime; 
  int nint   = ceil(fabs(t - tlo) * 2); // fabs is needed because we could
                                        // overshoot the root (by updating
                                        // g_starttime)
  int Ng     = 4;    // im using Ng = 4 because i have not implemented 5

  double vel = GaussQuad(Acceleration, tlo, t, nint, Ng);
  nfv += Ng*nint; 

  vel += g_startvel;

  // i think mignone wants this
  // (the number of function evaluations decreases a lot)
  g_startvel  = vel;
  g_starttime = t;

  cout << "tlo = " << tlo << "; thi = " << t << "; nint = " 
       << nint << "; Func eval = " << nfv << endl; 

  return vel;
}

double Acceleration(double t) {
  return (- exp(- t) + 1) / sqrt(t * t * t * t + 1);
} 
