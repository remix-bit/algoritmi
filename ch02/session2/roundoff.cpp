#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

int ex1(float, float&, float&, float&);
int ex2(float, float&, float&, float&);

int main() {

  cout << setiosflags(ios_base::scientific) << setprecision(7); 

  cout << "sqrt(x^2 + 1) - x" << endl;
  cout << "====================" << endl;
  
  float x;
  float fx1, fx2, ftayl; 

  for (int i = 0; i < 10; i++) {
    x = pow(10.0f, i+2);
    int flag = ex1(x, fx1, fx2, ftayl);

    cout << "x = " << x << "; fx1 = " << fx1 << "; fx2 = " << fx2 
         << "; f(taylor) = " << ftayl << endl;
  }

  cout << endl;
  cout << "1 - cos(x)" << endl;
  cout << "====================" << endl;


  for (int i = 0; i < 10; i++) {
    x = pow(10.0f, - i - 1);
    int flag = ex2(x, fx1, fx2, ftayl);

    cout << "x = " << x << "; fx1 = " << fx1 << "; fx2 = " << fx2 
         << "; f(taylor) = " << ftayl << endl;
  }

  return 0;
}

int ex1(float x, float& fx1, float& fx2, float& ftayl) {
  fx1 = sqrt(x * x + 1.0f) - x;
  fx2 = 1.0f / (sqrt(x*x + 1.0f) + x);

  // we need six significant digits... one or two terms are sufficient in order to reach
  // machine precision
  ftayl = 1 / (2.0f * x) - 1 / (8.0f * x * x * x);

  return 0;
}

int ex2(float x, float& fx1, float& fx2, float& ftayl) {
  fx1 = 1.0f - cos(x);
  fx2 = ( sin(x) * sin(x) ) / ( 1.0f + cos(x) );
  ftayl = ( x * x ) / ( 2.0f) - ( x * x * x * x ) / ( 24.0f );

  return 0;
}
