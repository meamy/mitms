#include "gate.h"
#include "ring.h"

int main() {
  Elt R(1, 2, 3, 4, 0);
  Elt A(0, 0, 0, 1, 1);
  Elt B = R+A;
  B.print();
  cout << "\n"; 
  test();
  return 0;
}
