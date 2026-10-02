#include <iostream>

using namespace std;

int main() {
  // it works but it could be better: when using non-defined values in
  // calculations the compiler may choose to cast everything as a double
  // now everything is best practice: we could also write 1.0f to enforce
  // interpretation as a float
  float x;
  float epsx = 1.0;
  float divFactor = 10.0;
  float x0 = 1.0;

  do {
    epsx /= divFactor;
    x = x0 + epsx; 
  } while ( x > x0 );

  // we only find the order of magnitude
  cout << "final eps for float = " << epsx << endl; 

  double y;
  double epsy = 1.0;
  double y0 = 1.0;

  do {
    epsy /= divFactor;
    y = y0 + epsy;
  } while ( y > y0 );

  cout << "final eps for double = " << epsy << endl;

  return 0;
}
