#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

double func(double x) {
  return sin(x);
}

double exder(double x) {
  return cos(x);
}

double ForwardDiff(double (*) (double), double, double);
double BackwardDiff(double (*) (double), double, double);
double CentralDiff(double (*) (double), double, double);
double FourthDiff(double (*) (double), double, double);

int main() {
  double h = 0.5;
  double x = 1;

  ofstream fdata;
  fdata.open("derivative.dat");

  double fd, bd, cd, fourthd, der;
  double errfd, errbd, errcd, errfourthd;

  for (int i = 1; i < 8; i++) {
    h *= 0.5;

    fd      = ForwardDiff(func, x, h);
    bd      = BackwardDiff(func, x, h);
    cd      = CentralDiff(func, x, h);
    fourthd = FourthDiff(func, x, h);

    der = exder(x); 

    errfd      = fabs(fd - der);
    errbd      = fabs(bd - der);
    errcd      = fabs(cd - der);
    errfourthd = fabs(fourthd - der);
    
    fdata << (double) 1 / h << "  " << errfd 
          << "  " << errbd << "  " << errcd <<  "  " << errfourthd << endl;
  }

  fdata.close();

  return 0;
}


double ForwardDiff(double (*func) (double), double x, double h) {
  return (func(x+h) - func(x)) / h;
}

double BackwardDiff(double (*func) (double), double x, double h) {
  return (func(x) - func(x-h)) / h;
}

double CentralDiff(double (*func) (double), double x, double h) {
  return (func(x+h) - func(x-h)) / (h * 2);
}

double FourthDiff(double (*func) (double), double x, double h) {
  return (func(x - h * 2) - func(x-h) * 8.0 + 
      func(x+h) * 8.0 - func(x + h * 2)) / (h * 12);
}
