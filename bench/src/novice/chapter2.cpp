#include <iostream>
#include "constants.h"

void compQuiz3() {
  using namespace std;
  cout << "Enter a double value: ";
  double x;
  cin >> x;

  cout << "Enter a second double value: ";
  double y;
  cin >> y;

  cout << "Enter one of the following: +, -, *, or / ";
  char op;
  cin >> op;

  if (op == '+')
    cout << x << " + " << y << " is " << x + y << endl;
  else if (op == '-')
    cout << x << " - " << y << " is " << x - y << endl;
  else if (op == '*')
    cout << x << " * " << y << " is " << x * y << endl;
  else if (op == '/')
    cout << x << " / " << y << " is " << x / y << endl;
}

double ballHeight(double init_h, double sec) {
  double fallen_h = constants::grav*sec*sec/2.0;
  double new_h = init_h-fallen_h;
  if (new_h < 0.0)
    new_h = 0.0;
  return new_h;
}

void compQuiz4() {
  using namespace std;
  cout << "Enter the height of a tower in meters: ";
  double x;
  cin >> x;

  cout << "At 0 seconds, the ball is at height: " << ballHeight(x,0) << endl;
  cout << "At 1 seconds, the ball is at height: " << ballHeight(x,1) << endl;
  cout << "At 2 seconds, the ball is at height: " << ballHeight(x,2) << endl;
  cout << "At 3 seconds, the ball is at height: " << ballHeight(x,3) << endl;
  cout << "At 4 seconds, the ball is at height: " << ballHeight(x,4) << endl;
  cout << "At 5 seconds, the ball is at height: " << ballHeight(x,5) << endl;
}

int main() {
  compQuiz3();
  compQuiz4();
  return 0;
}
