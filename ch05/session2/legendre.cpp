#include <iostream>
#include <iomanip>
#include <cmath>

#include "../../libs/root_finders.h"

#define DEBUG 0

using namespace std;

int g_legendredeg = 0; // = legendre polynomial/derivative degree
                       // func() and derfunc() take it implicitly
                       // in order to be double->double functions
                       // (it gets updated each time we change polynomial)

double legendrepol(double, int);
double derlegendrepol(double, int);

double func(double x) {
  return legendrepol(x, g_legendredeg);
}

double derfunc(double x) {
  return derlegendrepol(x, g_legendredeg);
}

int Bracket(double (*)(double), double, double, 
            int, double*, double*, int&);

int main() {
  double zero_bisection;
  double zero_falsepos;
  double zero_secant;
  double zero_newton;

  double bigA = -1;    // = left side of big interval
  double bigB = 1;     // = right side of big interval
  double xtol = 1e-7;
  double ytol = -1;

  cout << setiosflags(ios_base::scientific) << setprecision(12);

  // loop on the first 5 polynomial degrees
  for (g_legendredeg = 0; g_legendredeg < 6; g_legendredeg++) {
    cout << endl << "************************************" << endl;
    cout << "legendre with degree " << g_legendredeg << endl;

    const int n   = 21;                    // = init lenght of roots arrays
    int nbrackets = 2 * g_legendredeg + 1; // = number of brackets
                                           // (different from n for 
                                           // consistency with Mignone)
    int roots     = 0;
    double start[n];
    double end[n];

    Bracket(func, bigA, bigB, nbrackets,
                start, end, roots);
    cout << roots << " roots found" << endl << endl;

    double a, b; // = interval limits for each bracket
    for (int i = 0; i < roots; i++) {
      a = start[i];
      b = end[i];

      cout << endl << "interval [" << a << ", " << b << "]" << endl;

      Bisection(func, a, b, xtol, ytol, zero_bisection);
      FalsePosition(func, a, b, xtol, ytol, zero_falsepos);
      Secant(func, a, b, xtol, ytol, zero_secant);
      Newton(func, derfunc, a, b, xtol, ytol, zero_newton);

      cout << "bisection root      = " << zero_bisection << endl;
      cout << "false position root = " << zero_falsepos << endl;
      cout << "secant root         = " << zero_secant << endl;
      cout << "Newton root         = " << zero_newton << endl << endl;

      // computing weight for gauss integral algorithm
      // using the zeroes of the legendre polynomials
      double wi = (double) 2 / ((1.0 - zero_newton) 
                  * derfunc(zero_newton) * derfunc(zero_newton));
      cout << endl << "weight = " << wi << endl;
    }
    cout << "************************************" << endl;
  }
  return 0;
}

double legendrepol(double x, int n) {
  double p0 = 1;
  if (n == 0) return p0;

  double p1 = x;
  if (n == 1) return p1;

  double newp;
  for (int i = 2; i <= n; i++) {
    newp =  (((double) i * 2.0 - 1.0) * x * p1 
              - ((double) i - 1.0) * p0) / (double) i;
    p0 = p1;
    p1 = newp;
  }

  return newp;
}

double derlegendrepol(double x, int n) { // derivative
  double p0 = 1;
  if (n == 0) return 0;

  double p1 = x;
  if (n == 1) return 1;

  double newp;
  for (int i = 2; i <= n; i++) {        // same loop as legendre pol
    newp =  (((double) i * 2.0 - 1.0) * x * p1 
              - ((double) i - 1.0) * p0) / (double) i;
    p0 = p1;
    p1 = newp;
  }

  // given the last two polynomials p_n and p_(n-1) we can compute
  // the derivative of p_n
  return (x* p1 - p0) * (double) n / (x*x - 1.0); 
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
