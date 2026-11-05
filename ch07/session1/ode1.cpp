#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

using namespace std;

void EulerStep(double, double*, void (*) (double, double *, double *),
               double, int);
void dYdt(double, double*, double*);

double trueY(double t) {
  return exp(- t * t * 0.5);
}

int main() {
  const int neq = 1;

  double ts   = 0;
  double te   = 3;
  double Y[1] = {1};   // boundary value at ts

  int N    = 6;
  double h = te / N;
  double t;

  double correctsol;
  double abs_error, rel_error;
  cout << setiosflags(ios_base::scientific) << setprecision(6);
  
  // its better to extract a printout statement
  // than performing one step more?
  for (int i = 0; i < N+1; i++) {
    t = ts + h * i;
    correctsol = trueY(t);
    abs_error = fabs(Y[0] - correctsol);
    rel_error = abs_error / correctsol;
    cout << t << "  " << Y[0] << "  " << abs_error 
        << "  " << rel_error << endl;

    EulerStep (t, Y, dYdt, h, neq); 
  }

  return 0;
}


void EulerStep (double t, double *Y, 
    void (*RHSFunc)(double, double *, double *), double dt, int neq) {
  // Take one step dt using Euler method for the solution of dY/dt = rhs.  
  // Here neq is the number of ODE (the dimensionality of Y[]) 
  // and *RHSFunc() is a pointer to the actual function 
  // (in this case it points to dYdt()) that calculates the right-hand side
  // of the system of equations. 
  int k; 
  double rhs[256]; // Make sure rhs[] is large enough 
                   // (this implies neq < 256) 
                  // rhs[neq] is *NOT* standard C++ 
                  // (forbidden as static declaration) 
                  
  RHSFunc (t, Y, rhs); 
  for (k = 0; k < neq; k++) Y[k] += dt*rhs[k]; // Update solution array 
} 

// problem dependent (we are working with just one equation for now)
void dYdt (double t, double *Y, double *R) {
  double y = Y[0]; 
  R[0] = -t*y;
}
