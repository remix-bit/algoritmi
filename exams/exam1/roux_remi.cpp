// Name: Rémi Roux
// Date: 6 November 2026
//
// Code output:
// ******************************************************************
// Gauss L(n=10)  = 7.951983778230e+00
// Trapz L(n=10)  = 7.344182994898e+00
// Gauss L(n=100)  = 7.545534131724e+00
// Trapz L(n=100)  = 7.540217122615e+00
// Gauss L(n=1000)  = 7.540332706099e+00
// Trapz L(n=1000)  = 7.540331571809e+00
// ******************************************************************
#include <iostream>
#include <iomanip>
#include <cmath>

#define ECCENTRICITY 0.85
#define ALPHA 1
#define MAJORSEMIAXIS 1

using namespace std;

static double g_time = 0;

// root finder
int Secant(double (*)(double), double, double,
           double, double, double &);

// integrators
double TrapezoidQuad(double (*)(double), double, double, int);
double GaussQuad(double (*)(double), double, double, int, int); 

double vel(double);
double ecc_anomaly(double);
double lhs_ecc_anomaly(double);

int main() {
  // setting number output format
  cout << setiosflags(ios_base::scientific) << setprecision(12);

  int ng = 3;   // = gauss points
  int N  = 1;  //  = sub-interval points 

  double ts = 0;   // = start time
  double te = 10;  // = end time

  double trapezoid_integral, gauss_integral;

  for (int i = 0; i < 3; i++) {
    N *= 10;  // increasing number of sub-intervals for each iteration

    gauss_integral     = GaussQuad(vel, ts, te, N, ng);
    trapezoid_integral = TrapezoidQuad(vel, ts, te, N);

    cout << "Gauss L(n=" << N << ")  = " << trapezoid_integral << endl;
    cout << "Trapz L(n=" << N << ")  = " << gauss_integral << endl;
  }

  return 0;
}

// speed of planet at time t
double vel(double t) {
  // computing this factor to evaluate ecc_anomaly only once
  double factor = (ECCENTRICITY) * cos(ecc_anomaly(t));
  return (MAJORSEMIAXIS) * (ALPHA) * sqrt((factor + 1) / (- factor + 1));
}

// solving for eccentric anomaly with secant
double ecc_anomaly(double t) {
  double xtol = 1e-10;   // = tolerance
  double ytol = -1;
  double E;              // = root

  // graphically (but it can be shown precisely),
  // it is always true that:
  // t - ECCENTRICITY < E < t + ECCENTRICITY
  double a = t - 1;
  double b = t + 1;
  
  g_time = t;    // in order to use secant we need lhs function
                 // that depends on only one variable
  Secant(lhs_ecc_anomaly, a, b, xtol, ytol, E);

  return E;
}

// left hand side of equation for eccentric anomaly
// (the equation we want to solve is: lhs = 0)
double lhs_ecc_anomaly(double E) {
  return E - (ECCENTRICITY) * sin(E) - (ALPHA) * g_time;
}


// LIBRARY FUNCTIONS

int Secant(double (*func) (double), double a, double b, 
           double xtol, double ytol, double &zero) {
  double del    = b - a; // = x_k+1 - x_k (distance update)
  double olddel;         // = x_k - x_k-1 (needed to perform conv checks)
  int counter   = 1;

  // check if a is a root
  double fa = func(a);
  if (fa == 0) {
    zero = a;
    std::cout << "a = " << a << " is a zero: succesful secant" << std::endl;
    return 0;
  }

  // check if b is a root
  double fb = func(b);
  if (fb == 0) {
    zero = b;
    std::cout << "b = " << b << " is a zero: succesful secant" << std::endl;
    return 0;
  }

  while ((fabs(del) > xtol) && (fabs(fb) > ytol)) {
    if (counter > 128) {
      std::cout << "too many iterations: unsuccesful secant" << std::endl;
      return 1;
    }

    if (fabs(fb - fa) < 1e-18) {
      std::cout << "secant line is almost horizontal (div/0): "
           << "unsuccesful secant" << std::endl;
      return 1;
    }

    // applying secant algorithm
    counter++;
    olddel = del;                       // saving previous del
    del    = fb * (b - a) / (fb - fa);  // compute new del with rule
    a      = b;                         // a    <- x_k
    fa     = fb;                        // f(a) <- f(x_k)
    b     -= del;                       // b    <- x_k+1 = x_k - del
    fb     = func(b);                   // f(b) <- f(x_k+1)

    if (fabs(del) > fabs(olddel) * 5.0) {
      std::cout << "intervals started increasing in size: "
           << "unsuccesful secant" << std::endl;
      return 1;
    }
  }
  // saving last x_k+1 to zero
  zero = b;

  return 0;
}

double TrapezoidQuad(double (*F)(double), double a, double b, int N)  {
  double sum = 0.0;

  double h;
  h = (b - a) / (double)N;   // = interval width

  for (int i = 1; i < N; i++) {
    sum += F(a + i * h);     // sum 0.5 * F(x_i) + 0.5 * F(x_i)
                             // where x_i = a + i * h
                             // (every x_i delimits two intervals)
  }

  // sum 0.5 * F(a) + 0.5 F(b) 
  // (only point that delimit a single interval)
  sum += 0.5 * (F(a) + F(b));

  return sum * h;            // sum 0.5 * (F(x_i) + F(x_i+1)) * h
}

double GaussQuad(double (*F)(double), double a, double b, int N, int Ng) {
  double w[4], t[4]; // init of weigths and sample points for Ng < 5

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

    case 4: // the order of the pairs (point, param) does not matter
      t[0] = - sqrt((3. / 7.) - (2. / 7.) * sqrt(6. / 5.)) ;
      w[0] = (18. + sqrt(30.)) / 36.;
      t[1] = - t[0];
      w[1] = w[0];

      t[2] = - sqrt((3. / 7.) + (2. / 7.) * sqrt(6. / 5.));
      w[2] = (18. - sqrt(30.)) / 36.;
      t[3] = - t[2];
      w[3] = w[2];
      break;

    default:
      std::cout << "Unsupported value of Ng: " << Ng << std::endl;
      return 1.;
  }

  double sum  = 0.0;  // = integral across all intervals
  double sumk = 0.0;  // = (gauss) integral over each interval

  double x0, x1;      // = extremes for each interval

  double h = (b - a) / (double) N; // = width of intervals

  for (int i = 0; i < N; i++) { // looping over intervals
    x0 = a + i * h;             // init extremes of intervals
    x1 = a + (i + 1) * h;

    sumk = 0.0;

    for (int k = 0; k < Ng; k++) {  // looping over subintervals
      // adding contribution from each gauss interval ( h == x1 - x0)
      sumk += w[k] * F(0.5 * h * t[k] + 0.5 * (x1 + x0));
    }
    // adding each gaussian integral to total integral 
    sum += sumk;
  }

  return sum * 0.5 * h; // weighting by 0.5 * h at the end for performance
}
