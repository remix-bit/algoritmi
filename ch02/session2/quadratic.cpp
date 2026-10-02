#include <iostream>
#include <cmath>

using namespace std;

// returns exit code
int standardSolver(double, double, double, double&, double&);
int improvedSolver(double, double, double, double&, double&);

int main() {
  double x1 = 1.e-12;
  double x2 = 1.e12;
  double a = 1;
  double b = - (x1 + x2);
  double c = x1 * x2;

  double sol1, sol2;
  int flag = standardSolver(a, b, c, sol1, sol2);

  if (flag) { return 1; }

  cout << " --- standard solver --- " << endl;
  cout << "x1 = " << x1 << "; sol1 = " << sol1 << endl;
  cout << "x2 = " << x2 << "; sol2 = " << sol2 << endl;

  flag = improvedSolver(a, b, c, sol1, sol2);

  if (flag) { return 1; }

  cout << " --- improved solver --- " << endl;
  cout << "x1 = " << x1 << "; sol1 = " << sol1 << endl;
  cout << "x2 = " << x2 << "; sol2 = " << sol2 << endl;

  return 0;
}

int standardSolver(double a, double b, double c, double& sol1, double& sol2) {
  float sqDelta = sqrt(pow(b,2.0) - 4.0 * a * c);

  if (a == 0.0) {
    cout << "oops... division by zero!" << endl;
    return 1;
  }

  sol1 = (- b + sqDelta) / (2.0 * a);
  sol2 = (- b - sqDelta) / (2.0 * a);

  return 0;
}

// we are avoiding cancellation
int improvedSolver(double a, double b, double c, double& sol1, double& sol2) {
  float sqDelta = sqrt(pow(b,2.0) - 4.0 * a * c);

  if (a == 0.0) {
    cout << "oops... division by zero!" << endl;
    return 1;
  }

  if ( b < 0 ) {
    sol1 = (2.0 * c) / (- b + sqDelta);
    sol2 = (-  b + sqDelta) / (2.0 * a);
  } else {
    sol1 = (- b - sqDelta) / (2.0 * a);
    sol2 = (2.0 * c) / (- b - sqDelta);
  }

  return 0;
}
