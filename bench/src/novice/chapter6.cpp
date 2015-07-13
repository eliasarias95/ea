#include <iostream>

void section2Question1() {
  float hi_temp[365] = {0};
}

enum Animal {
  CHICKEN,
  DOG,
  CAT,
  ELEPHANT,
  DUCK,
  SNAKE,
  TOTAL
};

void section2Question2() {
  int legs[TOTAL] = {2,4,4,4,2,0};
}

void section3Question1() {
  int anArray[9] = {4,6,7,3,8,2,1,9,5};
  for (int i=0; i<9; ++i)
    std::cout << anArray[i] << " ";
  std::cout << "\n";
}

void section3Question2() {
  using namespace std;
  int nNumber;
  while (nNumber < 1 || nNumber > 9) {
    cout << "Enter a number: ";
    cin >> nNumber;
  }

  int anArray[9] = { 4, 6, 7, 3, 8, 2, 1, 9, 5 };
  for (int iii=0; iii < 9; iii++)
    cout << anArray[iii] << " ";
  cout << endl;
  for (int jjj=0; jjj< 9; jjj++) {
    if (anArray[jjj] == nNumber) {
      cout << "The number " << nNumber << " has index " << jjj << endl;
      break; // since each # in the array is unique, 
             //no need to search rest of array
    }
  }
}

int main() {
  //section2Question1();
  //section2Question2();
  //section3Question1();
  section3Question2();
  return 0;
}
