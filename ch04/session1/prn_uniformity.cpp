#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <cstdlib>

using namespace std;

int main() {
  srand48(time(NULL));

  ofstream fdata;  // declaring output stream 
                         
  fdata.open("prn_uniformity.dat"); // open output file 

  for (int i = 0; i < 1000; i++){ 
    fdata << i << "  " << drand48() << endl;
  }

  fdata.close();  

  fdata.open("prn_momenta.dat");
  
  for (int N = 1; N < 100000; N*=2) {
    double mean = 0.0;
    double mean_sq = 0.0;
    double x;
   
    for (int j = 0; j < N; j++) {
      x = drand48();
      mean += x;
      mean_sq += x * x;
    }
    
    mean = (mean / (double) N) - 0.5;
    mean_sq = (mean_sq / (double) N ) - ((double) 1.0 / (double) 3.0);
    fdata << N << "  " << fabs(mean) << "  " << fabs(mean_sq) << endl;
  }

  return 0;
}
