#include <iostream>
#include <cstdlib>
#include <ctime>

using namespace std;

int main() {
  srand(time(NULL));

  int down = 1;
  int up = 100;

  int number = rand() % up + 1;
  int guess;
  int nguess = 1;

    do {
    cout << "n in [" << down << ", " << up << "]" << endl; 
    cout << "type in your guess #" << nguess++ << " > ";
    cin >> guess;

    if (guess < number) {
      down = guess;
    } else if (guess > number) {
      up = guess;
    }
  } while (guess != number);

  cout << "good job! the number was found in " << 
    nguess - 1 << " iterations!" << endl;

  return 0;
}
