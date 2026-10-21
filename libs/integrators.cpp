#include "integrators.h"
#include <cmath>
#include <iostream>

double RectangularQuad(double (*F)(double), double a, double b, int N)  {
  double sum = 0.0;

  double h;
  h = (b - a) / (double)N;   // = interval width

  for (int i = 0; i < N; i++) {
    sum += F(a + i * h);    // sum F(x_i) where x_i = a + i * h
  }

  return sum * h;           // sum F(x_i) * h (are of each rectangle)
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

double SimpsonQuad(double (*F)(double), double a, double b, int N)  {
  double sum = 0.0;

  // N must be even
  if (N % 2 != 0) {
    std::cout << "N = " << N << ", but it should be even" << std::endl;
    return 1.;
  }

  double h;
  h = (b - a) / (double)N;   // = interval width

  for (int i = 1; i < N; i++) {
    if (i % 2 == 1) {             // odd elements are weighted by 4
      sum += (4.) * F(a + i * h); // 4 * F(x_i)
    } else {                      // even elements by 2
      sum += (2.) * F(a + i * h); // 2 * F(x_i)
    }
  }

  // the first and last elements are weighted by 1
  sum += F(a) + F(b);

  return sum * h / 3.;  // (h / 3) * (F_0 + 2 * F_1 + 4 * F_2 + ...  + F_N)
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

double GaussQuad2D(double (*F)(double, double), 
    double xa, double xb, 
    double ya, double yb, 
    int N, int Ng) {

  double w[4], t[4];  // init of weights and sample points for Ng < 5
                    
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

  double dx = (xb - xa) / (double) N; // = step in x direction
  double dy = (yb - ya) / (double) N; // = step in y direction

  double sum = 0.0; // = total integral

  // loops to define squares on which gaussian integration is performed
  for(int j = 0; j < N; j++) {
    for(int i = 0; i < N; i++) {
      double xc = xa + (0.5 + i) * dx; // = coords of small square centre
      double yc = ya + (0.5 + j) * dy;

      // gaussian integration loops
      double sumy = 0.0;
      for(int jk = 0; jk < Ng; jk++){ // for each y, gaussian on x
        double sumx = 0.0;
        for(int ik = 0; ik < Ng; ik++) {
          sumx += w[ik] * F(xc + 0.5 * dx * t[ik], yc + 0.5 * dy * t[jk]); 
        }
        sumx *= 0.5 * dx;
        sumy += w[jk]* sumx; // this is the gauss weighting on y
      }
      sumy *= 0.5 * dy;
      sum += sumy;
    }
  }
  return sum;
}
