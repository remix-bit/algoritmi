#include <iostream>
#include <cmath>
#include <iomanip> // to set number of significant digits

using namespace std;

int main() {
  double x;

  double xend = M_PI;
  int N = 778;
  double dx = xend / N;

  int steps = 0;

  for ( x = 0.0; x < xend; x += dx ) { steps++; }
  
  cout << "x with double loop method = " << x << endl;

  x = 0.0;

  // this is obviously better: we risk taking more steps (and stuff)
  for ( int i = 0; i < N; i++ ) { x += dx; }

  cout << "x with int loop method = " << x << endl;

  cout << "(because double used " << steps << " steps, while int " << N << " steps" << endl;

  return 0;
}
