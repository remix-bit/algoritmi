#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <ctime>

using namespace std;

int main() {
  srand48(time(NULL));

  ofstream fdata;

  int N = 2;
  int Ninside;
  double tol = 1e-5;

  double x;
  double y;
  double rel_err, err, pi_mc;

  fdata.open("pi.dat");

  // for the second point of the problem: you shouldn't restart the
  // integration for each N. simply add a point until the error reaches
  // the tolerance required.
  // the implementation here is worse.
  do {
    Ninside = 0;
    N      *= 2; // if you want better performance add a constant!
    for (int i = 0; i < N; i++) {
      x = drand48() * 2.0 - 1.0;
      y = drand48() * 2.0 - 1.0;
      if (x*x + y*y < 1) Ninside++;
    }

    pi_mc = 4.0 * (double) Ninside / (double) N;
    rel_err = fabs(pi_mc / M_PI - 1.0);
    err = fabs(pi_mc - M_PI);

    // to compute the exact version of the error ... TODO

    cout << "pi = " << pi_mc << " with rel_err = " << rel_err 
      << " with N = " << N << endl;
    fdata << N << "  " << err << endl;
  } while (rel_err > tol);

  fdata.close();

  return 0;
}
