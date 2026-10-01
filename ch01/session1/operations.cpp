#include <iostream>

using namespace std;

// he added controls for division by 0, and returned 1

int main() {
  float x = 5.3;
  float y = 1.2;

  int i = 20;
  int j = 6;

  // perform some operations ...

  cout << "--- float operations ---" << endl;
  cout << "x + y = " << x + y << endl;
  cout << "x - y = " << x - y << endl;
  cout << "x * y = " << x * y << endl;
  if (y != 0.0) cout << "x / y = " << x / y << endl;
  cout << endl;

  cout << "---- int operations ----" << endl;
  cout << "i + j = " << i + j << endl;
  cout << "i - j = " << i - j << endl;
  cout << "i * j = " << i * j << endl;
  if (j != 0) {
    cout << "i / j = " << i / j <<  " with remainder " << i % j << endl;
  }
  return 0;
}
