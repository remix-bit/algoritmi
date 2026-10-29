// Name: Rémi Roux
// Date: 24/10/2026 (30mins)
// Code output:
// ********************************************************
// Bisection, results:
// s = 1; Root: 7.12337528e+00
// s = 2; Root: 2.81527537e+01
// s = 3; Root: 6.14337274e+01
// Secant, results:
// s = 1; Root: 7.12337527e+00
// s = 2; Root: 2.81527537e+01
// s = 3; Root: 6.14337274e+01
// ********************************************************

#include <iostream>
#include <iomanip>
#include <cmath>

#include "../../libs/root_finders.h"

using namespace std;

double equation(double E, int s) {
  double V1  = 250;
  double V2  = 80;
  double Esq = sqrt(E);
  return Esq + asin(Esq / sqrt(V1))
             + asin(Esq / sqrt(V2)) - M_PI * (double) s;
}

double func1(double E) {
  return equation(E, 1);
}

double func2(double E) {
  return equation(E, 2);
}

double func3(double E) {
  return equation(E, 3);
}

int Bracket(double (*)(double), double, double,
            int, double*, double*, int&);

int main() {
  //double a = 0;
  //double b = 80;
  double xtol = 1e-8;
  double ytol = -1;

  cout << setiosflags(ios_base::scientific) << setprecision(8);

  double zero_bisection, zero_secant;
  double bisection[3], secant[3];

  double start[3] = {5, 25, 60};
  double end[3]   = {10, 30, 65};

  double (*funcs[3])(double) = {func1, func2, func3};

  for (int i = 0; i < 3; i++) {
    Bisection(funcs[i], start[i], end[i], xtol, ytol, zero_bisection);
    bisection[i] = zero_bisection;

    Secant(funcs[i], start[i], end[i], xtol, ytol, zero_secant);
    secant[i] = zero_secant;
  }

  // printing results to screen
  // (ntry for rootfinders TODO)
  cout << "Bisection, results:" << endl;
  for (int i = 0; i < 3; i++) {
    cout << "s = " << i + 1 << "; Root: " << bisection[i] << endl;
  }

  cout << "Secant, results:" << endl;
  for (int i = 0; i < 3; i++) {
    cout << "s = " << i + 1 << "; Root: " << secant[i] << endl;
  }

  return 0;
}
