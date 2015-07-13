#include <iostream>

int readNumber() {
  using namespace std;
  cout << "Enter an integer value: ";
  int x;
  cin >> x;
  return x;
}

void writeAnswer(int x) {
  using namespace std;
  cout << "Your answer is: " << x << endl;
}
