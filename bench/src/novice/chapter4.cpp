#include <iostream>

void section1aQuestion1() { //x and y created here
  using namespace std;
  cout << "Enter an integer: ";
  int x;
  cin >> x;

  cout << "Enter a second integer, larger than the first: ";
  int y;
  cin >> y;

  if (x>y) { // temp created here
    int temp = x;
    x = y;
    y = temp;
  } // temp destroyed here
  cout << "The smaller variable is " << x << endl;
  cout << "The larger variable is " << y << endl;
} // x and y destroyed here

void section4bQuestion1() {
  using namespace std;
  cout << "Enter your full name: ";
  string name;
  getline(cin,name);

  cout << "Enter your age: ";
  int x;
  cin >> x;

  cout << "You've lived " << static_cast<float>(x)/name.length() << 
    " years for each letter in your name.\n";
}

struct Advertising {
  int ads_shown;
  float ads_clicked;//as percentage
  float avg_earned;//from each ad clicked
};

void printAdStats(Advertising ad) {
  using namespace std;
  cout << "The number of ads shown is " << ad.ads_shown << "\n";
  cout << "The percentage of ads clicked is " << ad.ads_clicked << "\n";
  cout << "The average amount earned from each clicked ad is " << 
    ad.avg_earned << "\n";
  cout << "In total, you earned $" << 
    ad.ads_shown*ad.ads_clicked*ad.avg_earned << " today.\n";
}

struct Fraction {
  int num;
  int denom;
};

void multiply(Fraction fr1, Fraction fr2) {
  float dec1 = static_cast<float>(fr1.num)/static_cast<float>(fr1.denom);
  float dec2 = static_cast<float>(fr2.num)/static_cast<float>(fr2.denom);
  std::cout << "The product of these fractions is: " << dec1*dec2 << "\n";
}

void section7Question1() {
  using namespace std;
  Advertising ad;
  cout << "Enter the number of ads shown: ";
  cin >> ad.ads_shown;
  cout << "Enter the percentage of ads clicked: ";
  cin >> ad.ads_clicked;
  cout << "Enter the average amount earned from each clicked ad: ";
  cin >> ad.avg_earned;
  printAdStats(ad);
}

void section7Question2() {
  using namespace std;
  Fraction fr1;
  cout << "Enter the numerator to a fraction: ";
  cin >> fr1.num;
  cout << "Now enter the denominator: ";
  cin >> fr1.denom;
  Fraction fr2;
  cout << "Enter the numerator to a second fraction: ";
  cin >> fr2.num;
  cout << "Now enter the denominator: ";
  cin >> fr2.denom;
  multiply(fr1,fr2);
}

enum MonsterType {
  OGRE,
  DRAGON,
  ORC,
  GIANT_SPIDER,
  SLIME
};

struct Monster {
  MonsterType mt;
  std::string name;
  int health;
};

void printMonster(Monster mon) {
  std::string type;
  if (mon.mt == OGRE)
    type = "Ogre";
  else if (mon.mt == DRAGON)
    type = "Dragon";
  else if (mon.mt == ORC)
    type = "Orc";
  else if (mon.mt == GIANT_SPIDER)
    type = "Giant Spider";
  else if (mon.mt == SLIME)
    type = "Slime";
  else
    type = "Unknown";

  std::cout << "This " << type << " is named " << mon.name << " and has " << 
    mon.health << " health.\n";
}

void compQuiz() {
  Monster torg = {OGRE,"Torg",145};
  Monster blurp = {SLIME,"Blurp",23};
  printMonster(torg);
  printMonster(blurp);
}

int main() {
  //section4bQuestion1();
  //section1aQuestion1();
  //section7Question1();
  //section7Question2();
  compQuiz();
  return 0;
}
