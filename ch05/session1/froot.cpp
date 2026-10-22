#include <iostream>
#include <iomanip>
#include <cmath>

#define DEBUG 1 

using namespace std;

double func(double x) {
  return sin(x);
}

int Bisection(double (*)(double), double, double, double, double &);
int FalsePosition(double (*)(double), double, double, double, double &);

int main() {

  double zero_bisection;
  double zero_falsepos;

  double a = -3;
  double b = 8;
  double tol = 1e-7;

  cout << setiosflags(ios_base::scientific) << setprecision(6);

  Bisection(func, a, b, tol, zero_bisection);
  FalsePosition(func, a, b, tol, zero_falsepos);

  cout << "bisection root      = " << zero_bisection << endl;
  cout << "false position root = " << zero_falsepos << endl;

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

    if (fa * fb > 0) {
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

    if (fa * fb > 0) {
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

    if (fmid == 0) {
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
