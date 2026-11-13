#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

#define OMEGA 1

using namespace std;

void EulerStep(double, double*, void (*) (double, double *, double *),
               double, int);
void RK2Step(double, double*, void (*) (double, double *, double *),
               double, int);
void RK4Step(double, double*, void (*) (double, double *, double *),
               double, int);
void VerletPosStep(double*, double*, void (*) (double *, double *),
               double, int);
void dYdt(double, double*, double*);
void Acceleration(double*, double*);

double Energy(double x, double v) {
  double omegax = (OMEGA) * x;
  return (v * v + omegax * omegax) * 0.5;
}


int main() {
  const int neq = 2;

  // initial conditions
  double xs = 1;
  double vs = 0;

  ofstream fdata;

  double Y[2] = {xs, vs}; 

  double T  = M_PI * 2.0 / (OMEGA);      // = harmonic oscillator period
  double dt = T * 0.02;     // = time step
  double t  = 0;            // = live time

  int N = 10000;            // = max number of steps
  int Noscillations = 10;   // = max number of oscillations

  fdata.open("harmonic.dat");

  fdata << setiosflags(ios_base::scientific) << setprecision(6);

  fdata << t << "  " << Y[0] << "  " << Y[1] 
        << "  " << Energy(Y[0], Y[1]) << endl;


  for (int i = 0; i < N; i++) {
    RK2Step(t, Y, dYdt, dt, neq); 
    t += dt;

    fdata << t << "  " << Y[0] << "  " << Y[1] 
          << "  " << Energy(Y[0], Y[1]) << endl;

    if (t / T > Noscillations) break;
  }

  fdata << endl << endl;
  // initial conditions
  double pos[1] = {xs};
  double vel[1] = {vs};
  
  t = 0;              // = start time
  // dt = 0.02 * T;   // = time step

  const int nvar = 1;  // = number of variables

  fdata << t << "  " << pos[0] << "  " << vel[0] 
        << Energy(pos[0], vel[0]) << endl;

  for (int i = 0; i < N; i++) {
    VerletPosStep(pos, vel, Acceleration, dt, nvar); 
    t += dt;

    fdata << t << "  " << pos[0] << "  " << vel[0] 
          << "  " << Energy(pos[0], vel[0]) << endl;

    if (t / T > Noscillations) break;
  }

  fdata.close();

  return 0;
}

// problem dependent (we are working with just one equation for now)
void dYdt (double t, double *Y, double *R) {
  double x  = Y[0]; 
  double vx = Y[1];

  R[0] = vx;
  R[1] = - (OMEGA) * (OMEGA) * x;
}

void Acceleration (double *pos, double *acc) {
  double x = pos[0];

  acc[0] = - (OMEGA) * (OMEGA) * x;
}

void VerletPosStep(double *pos, double *vel,
    void (*acceleration)(double *, double *), double dt, int nvar) {
  for (int i = 0; i < nvar; i++) pos[i] += dt * 0.5 * vel[i];

  double accs[256];
  acceleration(pos, accs);

  for (int i = 0; i < nvar; i++) {
    vel[i] += dt * accs[i];
    pos[i] += dt * 0.5 * vel[i];
  }
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
  for (k = 0; k < neq; k++) Y[k] += dt*rhs[k]; // update solution array 
} 
