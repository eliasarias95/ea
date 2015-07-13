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


int main() {
  //section3Question1();
  section3Question2();
  return 0;
}
