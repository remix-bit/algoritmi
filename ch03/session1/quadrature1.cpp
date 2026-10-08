#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

double rectangular(double(*)(double), double, double, int);
double trapezoid(double(*)(double), double, double, int);
double simpson(double(*)(double), double, double, int);

// function on which we will apply the integration rules
// we pass this function through a pointer to the integration rules
double fun(double x) {
  return exp(-x);
}

int main() {
  double a = 0.0;
  double b = 1.0;
  int N = 4;

  double rect = rectangular(fun, a, b, N);
  double trap = trapezoid(fun, a, b, N);
  double simp = simpson(fun, a, b, N);

  cout << setiosflags(ios_base::scientific) << setprecision(13);
  cout << "Exact:       " << 1. - exp(-1.) << endl;
  cout << "Rectangular: " << rect << " (N = " << N << ")" << endl;
  cout << "Trapezoid:   " << trap << " (N = " << N << ")" << endl;
  cout << "Simpson:     " << simp << " (N = " << N << ")" << endl;

  cout << endl;
  cout << "------- now we iterate until convergence --------" << endl;
  double newrect, newtrap, newsimp, error;
  double tol = 1.e-5;
  int d = 1;

  do {
    d *= 2;
    newrect = rectangular(fun, a, b, d * N);
    error = fabs(newrect - rect);
    rect = newrect;
  } while(error > tol);

  cout << "Rectangular: " << newrect << " (N = " << d * N << ")" << endl;

  d = 1;
  do {
    d *= 2;
    newtrap = trapezoid(fun, a, b, d * N);
    error = fabs(newtrap - trap);
    trap = newtrap;
  } while(error > tol);

  cout << "Trapezoid:   " << newtrap << " (N = " << d * N << ")" << endl;

  d = 1;
  do {
    d *= 2;
    newsimp = simpson(fun, a, b, d * N);
    error = fabs(newsimp - simp);
    simp = newsimp;
  } while(error > tol);

  cout << "Simpson:     " << newsimp << " (N = " << d * N << ")" << endl;

  return 0;
}

double rectangular(double (*F)(double), double a, double b, int N)  {   
  double sum = 0.0;

  double h;
  h = (b - a) / (double)N;

  // multiply by h outside of the loop 
  for (int i = 0; i < N; i++) {
    sum += F(a + i * h);
  }

  return sum * h;
}


double trapezoid(double (*F)(double), double a, double b, int N)  {   
  double sum = 0.0;

  double h;
  h = (b - a) / (double)N;

  for (int i = 1; i < N; i++) {
    sum += F(a + i * h);
  }

  // first and last elements: excluding them is the best way (no ifs!)
  sum += 0.5 * (F(a) + F(b));

  return sum * h;
}

double simpson(double (*F)(double), double a, double b, int N)  {   
  double sum = 0.0;
  
  // we should check the parity of N
  
  double h;
  h = (b - a) / (double)N;

  for (int i = 1; i < N; i++) {
    if (i % 2 == 1) { // odd elements
      sum += (4.) * F(a + i * h);
    } else {         // even elements
      sum += (2.) * F(a + i * h);
    }
  }
  // the best way is defining w = 4, and then w = 6 - w at the end
  // to alternate for each iteration between 4 and 6
  
  // first and last elements
  sum += F(a) + F(b);

  // everything at the end
  return sum * h / 3.;
}
