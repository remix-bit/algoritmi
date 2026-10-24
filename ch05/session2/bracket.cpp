#include <iostream>
#include <iomanip>
#include <cmath>

#include "../../libs/root_finders.h"

#define DEBUG 0

using namespace std;

double func(double x) {
  return sin(x) - (x*x*0.01 + x * 0.2 + 1.0 / 3);
}

double derfunc(double x) { // derivative
  return cos(x) - (x * 0.02 + 0.2);
}

int Bracket(double (*)(double), double, double, 
            int, double*, double*, int&);

int main() {
  double zero_bisection;
  double zero_falsepos;
  double zero_secant;
  double zero_newton;

  double a    = -10;
  double b    = 10; 
  double xtol = 1e-7;
  double ytol = -1;

  cout << setiosflags(ios_base::scientific) << setprecision(6);

  const int n = 20;
  int roots   = 0;
  double start[n];
  double end[n];
  Bracket(func, a, b, n, start, end, roots);

  for (int i = 0; i < roots; i++) {
    a = start[i];
    b = end[i];

    cout << "root search in (" << a << ", " << b << ")" << endl;

    Bisection(func, a, b, xtol, ytol, zero_bisection);
    FalsePosition(func, a, b, xtol, ytol, zero_falsepos);
    Secant(func, a, b, xtol, ytol, zero_secant);
    Newton(func, derfunc, a, b, xtol, ytol, zero_newton);

    cout << "bisection root      = " << zero_bisection << endl;
    cout << "false position root = " << zero_falsepos << endl;
    cout << "secant root         = " << zero_secant << endl;
    cout << "Newton root         = " << zero_newton << endl << endl;
  }

  return 0;
}

int Bracket(double (*func) (double), double a, double b,
            int n, double *start, double *end, int &possibleroots) {
  double len = b - a;
  double del = len / n;

  double c = a;
  double d = a + del;

  double fc = func(c);
  double fd = func(d);

  possibleroots = 0;

  for (int i = 0; i < n + 1; i++) {
    #if DEBUG
    cout << "(c,d) = (" << c << ", " << d << ");  "
         << "(fc, fd) = (" << fc << ", " << fd  << ");  " << endl;
    #endif

    if (fc * fd < 0) {
      start[possibleroots] = c;
      end[possibleroots]   = d;
      possibleroots++;
      #if DEBUG
      cout << " ^ added!" << endl;
      #endif
    }

    c = a + del * i;
    d = a + del * (i + 1);

    fc = fd;
    fd = func(d);
  }

  return 0;
}
