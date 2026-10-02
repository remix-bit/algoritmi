#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

// great method: it converges with a square! (we'll understand how with Newton's method)
int main() {
  double S, x, xnew, error; // be attentive: one should define variables when init
  double tol = 1e-8; // tolerance

  cout << "enter real number: ";
  cin >> S;

  cout << "enter guess for sqrt: ";
  cin >> x;

  cout << setiosflags(ios_base::scientific) << setprecision(7);

  int counter = 1;
  do {
   xnew = 0.5 * ( x + ( S / x ) );
   error = fabs(xnew - x);

   cout << "iteration " << counter++ << ": x = " << xnew << "; err = " << error << endl;
   x = xnew;
  } while (error > tol);

  cout << "heron's sqrt( " << S << " ) = " << x << endl;
  cout << "true sqrt = " << sqrt(S) << endl;

  return 0;
}
