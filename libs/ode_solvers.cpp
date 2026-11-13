#include "ode_solvers.h"

void VerletPosStep(double *pos, double *vel,
    void (*acceleration)(double *, double *), double dt, int nvar) {
  // starting half step for positions
  for (int i = 0; i < nvar; i++) pos[i] += dt * 0.5 * vel[i];

  // computing accelerations after half step (sort of RHS)
  double accs[256];
  acceleration(pos, accs);

  // update velocities with full step, and positions with another half step
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
  double rhs[256]; 
                  
  // compute RHS (= derivative) for this step
  RHSFunc (t, Y, rhs); 

  // update Y with += RHS * dt
  for (int k = 0; k < neq; k++) Y[k] += dt*rhs[k]; 
} 
