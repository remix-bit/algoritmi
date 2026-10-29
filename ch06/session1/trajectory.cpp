#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

double g_alpha = 10;

double pos(double t) {
  if (t == 0) return 0.0;

  return g_alpha * t * t - t*t*t*
          ((double) 1 - exp(- g_alpha * g_alpha / t)); 
}

int main() {
  int N    = 100;
  double h = g_alpha / N;
  double t = 0;

  ofstream fdata;
  fdata.open("trajectory.dat");

  double vel, acc;
    


  for (int i = 1; i < N; i++) {

  }

  fdata.close();

  return 0;
}

