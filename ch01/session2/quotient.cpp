#include <iostream>

using namespace std;

void integerDivision(int, int, int &, int &);

int main() {
  int a = 8;
  int b = 3;

  int q, r;
  integerDivision(a, b, q, r);

  return 0;
}

void integerDivision(int a, int b, int & q, int & r) {
  if (b == 0) {
    cout << "Division by 0 is impossible!" << endl;
    exit(1); // this is worse: error control should be made inside of main
             // return an int to use as a flag!
  }

  q = a / b;
  r = a % b;
  cout << a << " / " << b << " = " << q << " with remainder "
       << r << endl;
}
