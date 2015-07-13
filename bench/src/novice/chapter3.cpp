#include <iostream>

bool isEven(int x) {
  return x%2==0;
}

void section2Question2() {
  using namespace std;
  cout << "Enter an integer value: ";
  int x;
  cin >> x;
  if (isEven(x))
    cout << "The number you entered is even.\n";
  else
    cout << "The number you entered is odd.\n";
}

void section7Question6() {
  using namespace std;
  cout << "Enter a number between 0 and 255: ";
  int x;
  cin >> x;
  if (x>= 0 && x <=255) {
    int y, count;
    while (x!=0) {
      y  = x%2;
      x /= 2;
      cout << y;
      ++count;
      if (count%4==0) cout << " ";
    }
    cout << endl;
  } 
  else
    cout << " " << endl;
}

int main() {
  section2Question2();  
  section7Question6();  
  return 0;
}
