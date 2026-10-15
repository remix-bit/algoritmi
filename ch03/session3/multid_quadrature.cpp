#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

double GaussQuad2D(double(*)(double, double), 
    double, double, double, double, int, int);

double func2D(double x, double y) {
  return (x*x*x*x)*(y*y) + 2.0*(x*x)*(y*y) - (x*x)*y + 2;
}

double disk(double x, double y) {
  if (sqrt(x*x + y*y) <= 1) {
    return 1;
  } else {
    return 0;
  }
}

int main() {
  
  cout << setiosflags(ios_base::scientific) << setprecision(7);

  double xa = - 1;
  double xb = + 1;
  double ya = - 1;
  double yb = + 1;

  int N = 1;
  int Ng = 4;

  double int2d = GaussQuad2D(func2D, xa, xb, ya, yb, N, Ng);

  cout << "2D gauss [N = " << N << ", Ng = " << Ng 
    << "] integration = " << int2d << endl;

  int Npi = 0;
  double pi, err;
  double pinew = 0.0;
  double tol = 1e-5;
  do {
  pi = GaussQuad2D(disk, xa, xb, ya, yb, ++Npi, Ng);
  err = fabs(pi - pinew);
  cout << "error in pi at iteration N = " << Npi << " is " << err << endl;
  pinew = pi;
  } while (err > tol);

  cout << "pi = " << pi << " with N = " << Npi << " iterations" << endl;

  Npi = 0;
  do {
  pi = GaussQuad2D(disk, xa, xb, ya, yb, ++Npi, Ng);
  err = fabs(pi - M_PI);
  cout << "error in pi at iteration N = " << Npi << " is " << err << endl;
  } while(err > tol);

  // it does not work well because gaussian integration supposes
  // that a polynomial can approximate sufficiently well the function,
  // but here we're trying to integrate a step function... 
  cout << "pi = " << pi << " with N = " << Npi << " iterations" << endl;

  return 0;
}

double GaussQuad2D(double (*F)(double, double), 
    double xa, double xb, 
    double ya, double yb, 
    int N, int Ng) {

  double w[4], t[4]; // they should be inited correclty: no runtime init
                     // (mignone prefers it this way)

  // defining weights and parameters for each choice of Ng (< 5)
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

  double dx = (xb - xa) / (double) N;
  double dy = (yb - ya) / (double) N;

  double sum = 0.0;

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
        sumy += w[jk]* sumx; // this is the weighting for gauss int on y
      }
      sumy *= 0.5 * dy;
      sum += sumy;
    }
  }
  return sum;
}
