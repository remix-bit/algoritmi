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

int main() {
  const int neq = 2;

  ofstream fdata;

  double ts   = 0;
  double te   = 3;
  int N    = 4;
  double t;

  double errEuler, errRK2, errRK4;

  fdata.open("convstudy.dat");

  cout << setiosflags(ios_base::scientific) << setprecision(6);

  for (int j = 0; j < 10; j++) {
    double YEuler[2] = {1, 0};   // boundary values at ts
    double YRK2[2]   = {1, 0};  
    double YRK4[2]   = {1, 0};   

    double h = te / N;

    for (int i = 1; i < N+1; i++) {
      EulerStep(t, YEuler, dYdt, h, neq); 
      RK2Step(t, YRK2, dYdt, h, neq); 
      RK4Step(t, YRK4, dYdt, h, neq); 

      t = ts + h * i;
    }

    errEuler = fabs(YEuler[0] - cos(t));
    errRK2   = fabs(YRK2[0] - cos(t));
    errRK4   = fabs(YRK4[0] - cos(t));
    fdata << N << "  " << errEuler << "  " << errRK2 << "  " << errRK4 << endl;

    N *= 2;
  }

  fdata.close();

  return 0;
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

// problem dependent (we are working with just one equation for now)
void dYdt (double t, double *Y, double *R) {
  double x = Y[0]; 
  double y = Y[1];
  R[0] = y;
  R[1] = -x;
}
