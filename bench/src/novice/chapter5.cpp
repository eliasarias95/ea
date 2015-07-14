#include <iostream>

int calculate(int x, int y, char op) {
  using namespace std;
  switch(op) {
    case '+':
      return x+y;
    case '-':
      return x-y;
    case '*':
      return x*y;
    case '/':
      return x/y;
    default:
      cout << "Error, invalid operator.\n";
      return 0;
  }
}

void section3Question1() {
  using namespace std;
  cout << "Enter an integer: ";
  int x;
  cin >> x;

  cout << "Enter a second integer: ";
  int y;
  cin >> y;

  cout << "Enter an operator (+, -, *, /): ";
  char op;
  cin >> op;

  cout << x << " " << op << " " << y << " = " << calculate(x,y,op) << "\n";
}

enum Animal {
  PIG,
  CHICKEN,
  GOAT,
  CAT,
  DOG,
  OSTRICH,
  TOTAL
};

std::string getAnimalName(Animal a) {
  switch(a) {
    case PIG:
      return "pig";
    case CHICKEN:
      return "chicken";
    case GOAT:
      return "goat";
    case CAT:
      return "cat";
    case DOG:
      return "dog";
    case OSTRICH:
      return "ostrich";
    default:
      std::cout << "Not a valid animal.\n";
  }
}

void printNumberOfLegs(Animal a) {
  std::cout << "A " << getAnimalName(a) << " has ";
  switch(a) {
    case PIG:
    case GOAT:
    case CAT:
    case DOG:
      std::cout << "4 legs.\n";
      break;
    case CHICKEN:
    case OSTRICH:
      std::cout << "2 legs.\n";
      break;
    default:
      std::cout << "Not a valid animal.\n";
  }
}

void section3Question2() {
  Animal cat = CAT;
  Animal chicken = CHICKEN;
  printNumberOfLegs(cat);
  printNumberOfLegs(chicken);
}

void section5Question2() {
  using namespace std;
  const int n = 26;
  char abc[n];
  for (int i=0; i<n; ++i) {
    abc[i] = 97+i;
    cout << abc[i] << " has ASCII code " << static_cast<int>(abc[i]) << "\n";
  }
}

/* other implementation
   void section5Question2() {
   char chValue = 'a';
   while (chValue <= 'z') {
   std::cout << chValue << " " << static_cast<int>(chValue) << "\n";
   chValue++;
   }
   }
   */

void section5Question3() {
  int outer = 5;
  while (outer >= 1) {
    int inner = outer;
    while (inner >= 1)
      std::cout << inner-- << " ";
    // print a newline at the end of each row
    std::cout << "\n";
    --outer;
  }
}

void section7Question1() {
  for (int i=0; i<=20; i+=2)
    std::cout << i << " ";
  std::cout << "\n";
}

int sumTo(int value) {
  int sum = 0;
  for (int i=0; i<=value; ++i)
    sum +=i;
  return sum;
}

void section7Question2() {
  std::cout << sumTo(5) << "\n";
}

int main() {
  //section3Question1();
  //section3Question2();
  //section5Question2();
  //section5Question3();
  //section7Question1();
  section7Question2();
  return 0;
}
