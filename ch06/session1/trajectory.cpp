#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

double g_alpha = 10;

double pos(double t) {
  if (t == 0) return 0.0;    // we cannot really represent it as a taylor
                             // series at t = 0... otherwise:
                             // if t < tol: return [series]

  return g_alpha * t * t - t*t*t*
          ((double) 1 - exp(- g_alpha * g_alpha / t)); 
}

int main() {
  int N        = 1000;
  double h     = g_alpha / N;
  double start = 0;

  ofstream fdata;
  fdata.open("trajectory.dat");

  double pos0, pos1, pos2;
  double vel, acc;
    
  // compute separately t = 0 case (discontinuous derivatives)

  double t = start;
  pos0 = pos(t);
  pos1 = pos(t + h);
  pos2 = pos(t + h * 2);
  
  vel = (pos1 - pos0) / h;                    // forward difference
  acc =  (pos2 + pos0 - pos1 * 2) / (h * h);  // adapted 2nd derivative

  fdata << t << "  " << pos0 << "  " << vel << "  " << acc << endl;

  for (int i = 1; i < N; i++) {
    t = start + h * i;
    vel = (pos2 - pos0) / (h *2);             // central difference
    acc = (pos2 + pos0 - pos1 * 2) / (h * h); // std 2nd derivative

    fdata << t << "  " << pos1 << "  " << vel << "  " << acc << endl;

    // update pos_(i-1), pos_i, pos_(i+1)
    pos0 = pos1;
    pos1 = pos2;
    pos2 = pos(start + h * (2 + i));
  }

  fdata.close();

  return 0;
}

