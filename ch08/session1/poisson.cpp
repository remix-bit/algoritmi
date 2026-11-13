#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "../../libs/root_finders.h"
#include "../../libs/ode_solvers.h"

using namespace std;

void dYdt(double, double*, double*);
double Residual(double);

int main() {
  // ------------------------------------------------------------- //
  // we begin by studying the function phi we get for various
  // initial guesses of s (the starting condition of phi')
  // ------------------------------------------------------------- //

  int N = 1000;  // = number of points we integrate 
  const int neq = 2; // = number of equations

  double a  = 0;       // = start point
  double b  = 20;      // = end point
  double dx = b / N;   // = step size
  double x;            // = indep live variable

  // opening file and setting output format
  ofstream fdata;
  fdata.open("poisson.dat");
  fdata << setiosflags(ios_base::scientific) << setprecision(6);

  double phis = 0;        // = initial condition for phi
  double s    = 0;        // = initial guess for derivative phi'
                          //   (initial condition)
  double Y[2] = {phis, s};

  for (int j = 0; j < 6; j++) {
    x = a;
    s = 0.2 * j;
    Y[0] = phis;
    Y[1] = s;

    fdata << x << "  " << Y[0] << "  " << Y[1] << endl;

    for (int i = 0; i < N; i++) {
      RK4Step(x, Y, dYdt, dx, neq);
      x += dx;

      fdata << x << "  " << Y[0] << "  " << Y[1] << endl;
    }

  fdata << endl << endl;
  }
  
  // ------------------------------------------------------------- //
  // now we plot the residual for s in [0,5]: graphically
  // we want to determine the correct initial guess
  // ------------------------------------------------------------- //

  double res;
  int M     = 100; 
  double S  = 0;
  double dS = 5.0 / 100;

  for (int j = 0; j < M; j++) {
    res = Residual(S); 
    fdata << S << "  " << res << endl;
    S += dS;
  }
  fdata << endl << endl;
  
  // ------------------------------------------------------------- //
  // now we find the exact value of s by using the bisection method
  // and then plot everything to compare against the exact solution
  // ------------------------------------------------------------- //

  double A    = 0;
  double B    = 1;
  double xtol = 1e-10;
  double ytol = -1;
  double exactS;

  Bisection(Residual, A, B, xtol, ytol, exactS);

  cout << "exact S = " << exactS << endl;

  x    = 0;
  Y[0] = 0;
  Y[1] = exactS;

  for (int i = 0; i < N; i++) {
    RK4Step(x, Y, dYdt, dx, neq);
    x += dx;

    fdata << x << "  " << Y[0] << "  " << Y[1] << endl;
  }

  fdata.close();

  return 0;
}

// problem dependent equation we want to solve (RHS)
void dYdt(double r, double *Y, double *R) {
  double phi    = Y[0];
  double derphi = Y[1];

  R[0] = derphi;
  R[1] = - r * 0.5 * exp(-r);
}

double Residual(double s) {
  int N = 1000;      // = number of points we integrate 
  const int neq = 2; // = number of equations

  double x  = 0;       // = start point (and live indep var)
  double b  = 20;      // = end point
  double dx = b / N;   // = step size

  // initial conditions
  double Y[2] = {0, s};
  for (int i = 0; i < N; i++) {
    RK4Step(x, Y, dYdt, dx, neq);
    x += dx;
  }
 
  return Y[0] - 1.0;
}
