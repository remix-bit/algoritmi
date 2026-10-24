#include <iostream>
#include <iomanip>
#include <cmath>

#define DEBUG 1 

using namespace std;

double func(double x) {
  return exp(-x) - x;
}

double derfunc(double x) { // derivative
  return - exp(-x) - 1;
}


int Bisection(double (*)(double), double, double, double, double &);
int FalsePosition(double (*)(double), double, double, double, double &);
int Secant(double (*)(double), double, double, double, double &);
int Newton(double (*)(double), double (*)(double), double, double, double, double &);

int main() {

  double zero_bisection;
  double zero_falsepos;
  double zero_secant;
  double zero_newton;

  double a = -1;
  double b = 1;
  double tol = 1e-7;

  cout << setiosflags(ios_base::scientific) << setprecision(6);

  Bisection(func, a, b, tol, zero_bisection);
  FalsePosition(func, a, b, tol, zero_falsepos);
  Secant(func, a, b, tol, zero_secant);
  Newton(func, derfunc, a, b, tol, zero_newton);

  cout << "bisection root      = " << zero_bisection << endl;
  cout << "false position root = " << zero_falsepos << endl;
  cout << "secant root         = " << zero_secant << endl;
  cout << "Newton root         = " << zero_newton << endl;

  return 0;
}

int Bisection(double (*func) (double), double a, double b,
              double tol, double &zero) {
  double fmid;
  int counter = 0;

  double fa = func(a);
  
  if (fa == 0) {
    zero = a;
    cout << "a = " << a << " is a zero: succesful bisection" << endl;
    return 0;
  }

  double fb = func(b);

  if (fb == 0) {
    zero = b;
    cout << "b = " << b << " is a zero: succesful bisection" << endl;
    return 0;
  }

  while ((fabs(b-a) > tol)) {
    if (counter > 128) {
      cout << "too many iterations: unsuccesful bisection" << endl;
      return 1;
    }

    if (fa * fb > 0) { // maybe it is better to check this only once outside
      cout << "an interval with f(a) * f(b) > 0 was reached: " 
           << "unsuccesful bisection" << endl;
      return 1;
    } 

    zero = 0.5 * (a + b);
    fmid = func(zero);

    counter++;
    #if DEBUG
    cout << "Bisection(): k = " << counter
         << "; [a, b] = [" << a << ", " << b << "]; " 
         << "xm = " << zero
         << "; dx = " << fabs(b-a)
         << "; fm = " << fmid << endl;
    #endif

    if (fmid == 0) {
      return 0;
    } else if (fa * fmid < 0) {
      b = zero;
      fb = fmid;
    } else {  // fa * fmid > 0 is the only possibility
      a = zero;
      fa = fmid;
    }
  }

  return 0;
}

int FalsePosition(double (*func) (double), double a, double b,
              double tol, double &zero) {
  double del = fabs(b-a);
  double fmid;
  int counter = 0;

  double fa = func(a);
  
  if (fa == 0) {
    zero = a;
    cout << "a = " << a << " is a zero: succesful false position" << endl;
    return 0;
  }

  double fb = func(b);

  if (fb == 0) {
    zero = b;
    cout << "b = " << b << " is a zero: succesful false position" << endl;
    return 0;
  }

  do {
    if (counter > 128) {
      cout << "too many iterations: unsuccesful false position" << endl;
      return 1;
    }

    if (fa * fb > 0) { // maybe it's better to check only once outside
      cout << "an interval with f(a) * f(b) > 0 was reached: " 
           << "unsuccesful false position" << endl;
      return 1;
    } 

    zero = (a * fb - b * fa) / (fb - fa);
    fmid = func(zero);

    counter++;
    #if DEBUG
    cout << "FalsePosition(): k = " << counter
         << "; [a, b] = [" << a << ", " << b << "]; " 
         << "xm = " << zero
         << "; |del| = " << del
         << "; fm = " << fmid << endl;
    #endif

    if (fmid == 0) { // i dont know if it makes sense...
      return 0;
    } else if (fa * fmid < 0) {
      del = fabs(b - zero);
      b  = zero;
      fb = fmid;
    } else { // fa * fmid > 0 is the only possibility
      del = fabs(a - zero);
      a  = zero;
      fa = fmid;
    }
  } while (del > tol);

  return 0;
}

int Secant(double (*func) (double), double a, double b, 
           double tol, double &zero) {
  double del = b - a;
  double olddel;
  int counter = 0;

  double fa = func(a);
  
  if (fa == 0) {
    zero = a;
    cout << "a = " << a << " is a zero: succesful secant" << endl;
    return 0;
  }

  double fb = func(b);

  if (fb == 0) {
    zero = b;
    cout << "b = " << b << " is a zero: succesful secant" << endl;
    return 0;
  }

  while (fabs(del) > tol) {
     if (counter > 128) {
      cout << "too many iterations: unsuccesful secant" << endl;
      return 1;
    }

    olddel = del;
    del = fb * (b - a) / (fb - fa);

    if (fabs(del) > fabs(olddel)) {
      cout << "intervals started increasing in size: "
           << "unsuccesful secant" << endl;
      return 1;
    }

    counter++;
    #if DEBUG
    cout << "Secant(): k = " << counter
         << "; xa = " << a << ", xb = " << b 
         << "; dx = " << del << endl;
    #endif

    a  = b;
    fa = fb;
    b  = b - del;
    fb = func(b);
  }

  zero = b;

  return 0;
}

int Newton(double (*func) (double), double(*derivative) (double),
           double a, double b, 
           double tol, double &zero) {
  int counter = 0;

  zero = 0.5 * (a + b);
  double fxc = func(zero);

  double del = fxc / derivative(zero);
 
  while (fabs(del) > tol) {
     if (counter > 128) {
      cout << "too many iterations: unsuccesful Newton" << endl;
      return 1;
    }

    if (fxc == 0) {
      cout << "xc = " << zero << " is a zero: succesful Newton" << endl;
      return 0;
    }

    del   = fxc / derivative(zero);
    zero -= del;
    fxc   = func(zero);

    counter++;
    #if DEBUG
    cout << "Newton(): k = " << counter
         << "; xc = " << zero 
         << "; dx = " << del << endl;
    #endif
  }

  return 0;
}
