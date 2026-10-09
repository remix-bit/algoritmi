#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

double simpson(double(*)(double), double, double, int);
double trapezoid(double(*)(double), double, double, int);
double gauss(double(*)(double), double, double, int, int);

double Si(double x) {
  return sin(x) / x;
}

double taylorSi(double x) {
  if (fabs(x) < 1.e-6) { // ci accolliamo questo if dai
    return 1. - (x*x) / 6. + (x*x*x*x) / 120.;
  } else {
    return sin(x) / x;
  }
}

int main() {
  double a = 0.0;
  double b = 0.8;
  int N = 1;
  int Ng = 3;

  cout << setiosflags(ios_base::scientific) << setprecision(7);

  cout << "[Gaussian n = " << N << ", ng = " << Ng << "] = " << 
    gauss(Si, a, b, N, Ng) << endl;

  for (int n = 1; n < 9; n *= 2) {
    cout << "[Trapezoidal n = " << n << "]      = " << 
      trapezoid(taylorSi, a, b, n) << endl;
  }

  for (int n = 2; n < 9; n *= 2) {
    cout << "[Simpson n = " << n << "]          = " << 
      simpson(taylorSi, a, b, n) << endl;
  }

  return 0;
}

double gauss(double (*F)(double), double a, double b, int N, int Ng) {
  double w[8], t[8]; // they should be inited correclty: no runtime init
                     // (mignone prefers it this way)

  // defining weights and parameters for each choice of Ng (< 9)
  switch (Ng) {
    case 1:
      w[0] = 2.0; t[0] = 0.0;
      break;

    case 2:
      w[0] = 1.0; t[0] = - sqrt(1./3.);
      w[1] = 1.0; t[1] = - t[0];
      break;

    case 3:
      w[0] = 5./9.; t[0] = - sqrt(0.6);
      w[1] = 8./9.; t[1] = 0.0;
      w[2] = w[0];  t[2] = - t[0];
      break;

    case 4: // it may be wrong check it AAAAAAAA
      w[0] = - sqrt((3. / 7.) - (2. / 7.) * sqrt(6. / 5.)) ; 
      t[0] = (18. + sqrt(30.)) / 36.;
      w[1] = - w[0];
      t[1] = t[0];

      w[2] = - sqrt((3. / 7.) + (2. / 7.) * sqrt(6. / 5.)); 
      t[2] = (18. - sqrt(30.)) / 36.;
      w[3] = - w[2];
      t[3] = t[2];

    default:
      cout << "Unsupported value of Ng: " << Ng << endl;
      return 1.;
}

  double sum  = 0.0;
  double sumk = 0.0;

  double x0, x1;     // extremes for each interval 
  
  double h = (b - a) / (double) N; // defining width of intervals

  for (int i = 0; i < N; i++) { // looping over total intervals
    x0 = a + i * h;
    x1 = a + (i + 1) * h;

    sumk = 0.0;

    for (int k = 0; k < Ng; k++) { 
      // adding contribution from each gauss interval ( h == x1 - x0)
      sumk += w[k] * F(0.5 * h * t[k] + 0.5 * (x1 + x0));
    }
    // adding to total integral (weighted at the end)
    sum += sumk;
  }

  return sum * 0.5 * h;
}

double simpson(double (*F)(double), double a, double b, int N)  {   
  double sum = 0.0;
  
  // we should check the parity of N
  if (N % 2 != 0) {
    cout << "N = " << N << ", but it should be even" << endl;
    return 1.;
  }
  
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
