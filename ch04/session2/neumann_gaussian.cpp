#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <cmath>

using namespace std;
  
double gauss(double x, double mu, double sigma) {
  return (1.0) / (sqrt(2.0 * M_PI) * 0.5) * 
    exp(-0.5 * (x - mu) * (x - mu) / (sigma * sigma));
}

int main() {
  srand48(time(NULL));
  ofstream fdata;

  int N = 1e5;
  double mu = 0.0;
  double sigma = 0.5;

  double C = 1 / sigma;
  double xmin = - 5.0;
  double xmax = + 5.0;
  double deltax = xmax - xmin;

  double x, y;
  fdata.open("neumann_gaussian.dat");
  for (int i = 0; i < N; i++) {
     x = xmin + drand48() * deltax;
     y = drand48() * C;

     if (y < gauss(x, mu, sigma)) {
      fdata << x << "  " << y << endl;
      }
  }

  fdata.close();

  return 0;
}

