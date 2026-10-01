#include <iostream>

using namespace std;

void printWhile() {
  int i = 1;
  while (i < 11) {
    cout << i++ << endl;
    // if (i % 2 == 1) cout << i << endl;
    // i++
  }
}

void printFor() {
  for (int i = 1; i < 11; i++) {
    cout << i << endl;
    // if (i % 2 == 1) cout << i << endl;
  }
}

int main() {
  cout << "printing 1 to 10 with while loop..." << endl;
  printWhile();
  cout << endl;

  cout << "printing 1 to 10 with for loop..." << endl;
  printFor();

  return 0;
}
