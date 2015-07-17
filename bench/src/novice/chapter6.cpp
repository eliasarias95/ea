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

/**
 * A pointer (pnPtr) can be set equal to an array (i.e. point to the 
 * beginning of the array) (e.g. pnPtr = szName). By adding integers to the
 * array, you are accessing different elements of the array 
 * (e.g. pnPtr = szName + 1 is equivalent to szName[1]).
 */
void section8() {
  const int nArraySize = 7;
  char szName[nArraySize] = "Mollie";
  int nVowels = 0;
  for (char *pnPtr = szName; pnPtr < szName + nArraySize; ++pnPtr) {
    switch (*pnPtr) {
      case 'A':
      case 'a':
      case 'E':
      case 'e':
      case 'I':
      case 'i':
      case 'O':
      case 'o':
      case 'U':
      case 'u':
        ++nVowels;
        break;
    }
  }
  std::cout << "szName+nArraySize = " << (szName) << "\n";
  std::cout << szName << " has " << nVowels << " vowels\n";
}

int main() {
  //section2Question1();
  //section2Question2();
  //section3Question1();
  //section3Question2();
  section8();
  return 0;
}
