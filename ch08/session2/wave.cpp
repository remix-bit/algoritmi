#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "../../libs/root_finders.h"
#include "../../libs/ode_solvers.h"

using namespace std;

double g_k;  // global definition of k 

void dYdt(double, double*, double*);

int main() {
  double xs = 0;   // = start point
  double xe = 1;   // = end point

       g_k = 1;    // = initial guess for eigenvalue
  double s = 1;    // = guess for first derivative 
                   //   (does not matter)

  int N     = 100;
  double dx = fabs(xe - xs) / N;
  const int neq = 2;

  ofstream fdata;
  fdata.open("wave.dat");
  fdata << setiosflags(ios_base::scientific) << setprecision(6);

  double Y[2] = {xs, s};
  double x = xs;

  fdata << x << "  " << Y[0] << "  " << Y[1] << endl;

  for (int i = 0; i < N; i++) {
    RK4Step(x, Y, dYdt, dx, neq);
    x += dx;
    fdata << x << "  " << Y[0] << "  " << Y[1] << endl;
  }

  fdata.close();

  return 0;
}

// problem dependent equation we want to solve (RHS)
void dYdt(double x, double *Y, double *R) {
  double phi    = Y[0];
  double derphi = Y[1];

  R[0] = derphi;
  R[1] = - g_k * g_k * phi;
}

