/* assert example */
#include <stdio.h>      /* printf */
#include <assert.h>     /* assert */
#include <iostream>

void print_number(int* myInt) {
  assert (myInt!=NULL && *myInt < 5);
  printf ("%d\n",*myInt);
}

void trace(std::string s) {
  std::cout << s << "\n";
}

int main ()
{
  int a=10;
  int * b = NULL;
  int * c = NULL;

  b=&a;

  //print_number (b);
  //print_number (c);
  trace("a= "+a);

  return 0;
}
