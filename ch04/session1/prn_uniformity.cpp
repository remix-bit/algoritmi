#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <cstdlib>

using namespace std;

int main() {
  srand48(time(NULL));

  ofstream fdata;        // declare Output stream class to operate on files 
                         
  fdata.open("prn_uniformity.dat"); // open output file 

  for (int i = 0; i < 1000; i++){ 
    fdata << i << "  " << drand48() << endl;
  }
  fdata.close();  

  int N = 10;
  fdata.open("prn_moments.dat");
  for (int i = 1; i < 18; i++) {
    double mean = 0.0;
    double mean_sq = 0.0;
    double x;
    for (int j = 0; j < N; j++) {
      x = drand48();
      mean += x;
      mean_sq += x * x;
    }
    mean = (mean / (double) N) - 0.5;
    mean_sq = (mean_sq / (double) N ) - (double) 1.0 / (double) 3.0;
    fdata << N << "  " << mean << "  " << mean_sq << endl;
    N *= 2;
  }

  return 0;
}
