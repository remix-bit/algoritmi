
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <cstdlib>

using namespace std;

int main() {
  srand48(time(NULL));

  ofstream fdata;

  fdata.open("decay.dat");

  int N = 1000; // = total number of particles
  int T = 900;  // = total time
  int dt = 1;   // = time step

  double lambda = 0.01;

  for (int t = 0; t < T; t += dt) {
    int decay = 0;
    for (int j = 0; j < N; j++) {
      if (drand48() < lambda) decay++;
    }
    N -= decay;

    fdata << t << "  " << N << endl;

    if (N == 0) break; // best practice put <?
  }

  fdata.close();

  return 0;
}
