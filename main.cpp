#include "gate.h"
#include "ring.h"

int main() {
  Ring R(1, 2, 3, 4, 0);
  Ring A(0, 0, 0, 1, 1);
  Ring B = R+A;
  B.print();
  cout << "\n"; 
  test();
  return 0;
}
