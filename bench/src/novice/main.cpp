/* assert example */
#include <stdio.h>      /* printf */
#include <assert.h>     /* assert */
#include <ucsl.h>
#include <iostream>

void print_number(int* myInt) {
  assert (myInt!=NULL && *myInt < 5);
  printf ("%d\n",*myInt);
}

int main ()
{
  int n  = 10;
  int *b = mem_alloc(n,sizeof(int));
  int *c = mem_alloc(n,sizeof(int));

  memset(b,1,n*sizeof(int));
  c = b;
  std::cout << c[4] << "\n";
  //print_number (b);
  //print_number (c);

  return 0;
}
