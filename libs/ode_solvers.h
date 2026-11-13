void EulerStep(double, double*, void (*) (double, double *, double *),
               double, int);
void RK2Step(double, double*, void (*) (double, double *, double *),
               double, int);
void RK4Step(double, double*, void (*) (double, double *, double *),
               double, int);
void VerletPosStep(double*, double*, void (*) (double *, double *),
               double, int);
