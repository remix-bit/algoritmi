#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

using namespace std;

void EulerStep(double, double*, void (*) (double, double *, double *),
               double, int);
void RK2Step(double, double*, void (*) (double, double *, double *),
               double, int);
void RK4Step(double, double*, void (*) (double, double *, double *),
               double, int);
void dYdt(double, double*, double*);

double orbitvel(double x, double y, double alpha){
  if (alpha >= 2) {
    cout << "attention! unbounded orbits for a = " << alpha << endl;
  }
  double r = sqrt(x * x + y * y);
  return sqrt(alpha / r);
}

double stepsize(double dtheta, double * Y) {
  double x  = Y[0];
  double y  = Y[1];
  double vx = Y[2];
  double vy = Y[3];
  return dtheta * sqrt((x * x + y * y) / (vx * vx + vy * vy));
}

int main() {
  const int neq = 4;

  // initial conditions
  double ecc = 0.3;
  double xs  = 4;
  double ys  = 0;
  double vxs = 0;
  double vys = orbitvel(xs, ys, ecc); 

  ofstream fdata;

  double Y[4] = {xs, ys, vxs, vys}; 

  double dtheta = 0.1;     // = angle step
  double h;                // = time step (inited later) 
  double ts     = 0;       // = start time
  double t      = ts;      // = live time

  int counterturns = 0;    // = counter for half-orbits
  int Norbits      = 10;   // = number of orbits in simulation
  double oldvx;
  int N = 100000;          // = max number of steps

  fdata.open("kepler.dat");

  fdata << setiosflags(ios_base::scientific) << setprecision(6);

  fdata << Y[0] << "  " << Y[1] << "  " << Y[2] << "  " << Y[3] << endl;

  for (int i = 0; i < N; i++) {
    h = stepsize(dtheta, Y); // defining time step for const dtheta

    oldvx = Y[2];
    RK4Step(t, Y, dYdt, h, neq); 

    t = ts + h * i;
    fdata << Y[0] << "  " << Y[1] << "  " << Y[2] << "  " << Y[3] << endl;

    if (oldvx * Y[2] <= 0) counterturns++;
    if (counterturns > 2 * Norbits) break;
  }

  cout << "turns = " << counterturns << endl;

  fdata.close();

  return 0;
}

// problem dependent (we are working with just one equation for now)
void dYdt (double t, double *Y, double *R) {
  double x  = Y[0]; 
  double y  = Y[1];
  double vx = Y[2];
  double vy = Y[3];

  double rcubed = sqrt(x*x + y*y);
  rcubed = rcubed * rcubed * rcubed;

  R[0] = vx;
  R[1] = vy;
  R[2] = - x / rcubed;
  R[3] = - y / rcubed;
}

void RK2Step (double t, double *Y, 
    void (*RHSFunc)(double, double *, double *), double dt, int neq) {
  double Y1[256], k1[256], k2[256];                

  RHSFunc(t, Y, k1); 
  for (int i = 0; i < neq; i++) Y1[i] = Y[i] + 0.5 * dt * k1[i]; 

  RHSFunc(t + 0.5 * dt, Y1, k2);
  for (int i = 0; i < neq; i++) Y[i] += dt * k2[i];
}

void RK4Step (double t, double *Y, 
    void (*RHSFunc)(double, double *, double *), double dt, int neq) {
  double Y1[256], Y2[256], Y3[256];
  double k1[256], k2[256], k3[256], k4[256];

  RHSFunc(t, Y, k1);
  for (int i = 0; i < neq; i++) Y1[i] = Y[i] + 0.5 * dt * k1[i];

  RHSFunc(t + 0.5 * dt, Y1, k2);
  for (int i = 0; i < neq; i++) Y2[i] = Y[i] + 0.5 * dt * k2[i];

  RHSFunc(t + 0.5 * dt, Y2, k3);
  for (int i = 0; i < neq; i++) Y3[i] = Y[i] + dt * k3[i];

  RHSFunc(t + dt, Y3, k4);
  for (int i = 0; i < neq; i++) {
    Y[i] += dt * (k1[i] + k2[i] * 2 + k3[i] * 2 + k4[i]) / 6;
  }
}

void EulerStep (double t, double *Y, 
    void (*RHSFunc)(double, double *, double *), double dt, int neq) {
  int k; 
  double rhs[256]; // Make sure rhs[] is large enough 
                   // (this implies neq < 256) 
                   // rhs[neq] is *NOT* standard C++ 
                   // (forbidden as static declaration) 
                  
  RHSFunc (t, Y, rhs); 
  for (k = 0; k < neq; k++) Y[k] += dt*rhs[k]; // Update solution array 
} 


