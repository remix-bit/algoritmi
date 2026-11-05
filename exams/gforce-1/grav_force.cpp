// Name: Rémi Roux
// Date: 5 November 2026 (45mins)
//
// Code output:
// ******************************************************************
//
// ******************************************************************
#include <iostream>
#include <iomanip>
#include <cmath>
#include "../../libs/integrators.h"

using namespace std;

double g_b = M_PI;

double pos(double);
double integrand_x(double);
double integrand_y(double);

int main() {
  double Wx, Wy; // work in x and y directions
  double W    = 0; 
  double oldW = 0;        

  double tol = 1e-6;
  double err = 1;

  double nsub = 2;

  double a = 0;
  double b = M_PI * 4;

  cout << setiosflags(ios_base::scientific) << setprecision(6);

  do {
    nsub *= 2;
    Wx = TrapezoidQuad(integrand_x, a, b, nsub);
    Wy = TrapezoidQuad(integrand_y, a, b, nsub);
    oldW = W;
    W = Wx + Wy;
    err = fabs(W - oldW);
  } while (err > tol);

  cout << "Trapezoidal: n = " << nsub << "; Wx = " << Wx << "; Wy = " 
       << Wy << "; Wtot = Wx + Wy = " << W << endl;

  nsub      = 2;
  double ng = 3;
  W         = 0;
  oldW      = 0;

  do {
    nsub *= 2;
    Wx = GaussQuad(integrand_x, a, b, nsub, ng);
    Wy = GaussQuad(integrand_y, a, b, nsub, ng);
    oldW = W;
    W = Wx + Wy;
    err = fabs(W - oldW);
  } while (err > tol);

  cout << "Gauss:       n = " << nsub << "; Wx = " << Wx << "; Wy = " 
       << Wy << "; Wtot = Wx + Wy = " << W << endl;

  // now i work out the exact value

  double ri = pos(a);
  double rf = pos(b);
  double trueW = - (1.0 / ri) + (1.0 / rf);   // - Delta U

  cout << "(Exact = " << trueW << ")" << endl;

  return 0;
}

double pos(double theta) {
  double toverb = theta / g_b;
  return (toverb + 1) * exp(- toverb);
}

double integrand_x(double theta) {
  double r = pos(theta);
  double toverb = theta / g_b;
  double exptoverb = exp(- toverb);

  double drdt = - toverb * exptoverb / g_b;
  double dxdt = drdt * cos(theta) - r * sin(theta);

  double integrand = (- cos(theta) / (r * r)) * dxdt;

  return integrand;
}

double integrand_y(double theta) {
  double r = pos(theta);
  double toverb = theta / g_b;
  double exptoverb = exp(- toverb);

  double drdt = - toverb * exptoverb / g_b;
  double dydt = drdt * sin(theta) + r * cos(theta);

  double integrand = (- sin(theta) / (r * r)) * dydt;

  return integrand;
}
