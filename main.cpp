#include "gate.h"
#include "matrix.h"

int main() {
  Elt R(1, 2, 3, 4, 0);
  Elt A(0, 0, 0, 1, 1);
  Elt B = R.conj();
  Rmatrix M(2, 2);
  M(0, 0) = Elt(1, 0, 0, 0, 0);
  M(0, 1) = Elt(0, 0, 0, 0, 0);
  M(1, 0) = Elt(0, 0, 0, 0, 0);
  M(1, 1) = Elt(1, 0, 0, 0, 0);
  M.print();
  cout << "\n"; 
  test();
  return 0;
}
