#include "ring.h"

#define PI 3.14159

Elt::Elt() {a = b = c = d = n = 0;}
Elt::Elt(int aa, int bb, int cc, int dd, int nn) { 
  a = aa; b = bb; c = cc; d = dd; n = nn;
}
Elt::Elt(const Elt & R) {
  *this = R;
}

Elt & Elt::operator=  (const Elt & R) {
  a = R.a; b = R.b; c = R.c; d = R.d; n = R.n;
}
Elt & Elt::operator+= (const Elt & R) {
  const Elt * low, * high;
  int diff;
  if (n < R.n) {
    low = this;
    high = & R;
    diff = pow(2, R.n - n);
  } else {
    low = & R;
    high = this;
    diff = pow(2, n - R.n);
  }

  a = diff*low->a + high->a;
  b = diff*low->b + high->b;
  c = diff*low->c + high->c;
  d = diff*low->d + high->d;
  n = high->n;
}
Elt & Elt::operator-= (const Elt & R) {
  const Elt * low, * high;
  int s1, s2;
  if (n < R.n) {
    low = this;
    high = & R;
    s1 = pow(2, R.n - n);
    s2 = -1;
  } else {
    low = & R;
    high = this;
    s1 = -pow(2, n - R.n);
    s2 = 1;
  }

  a = s1*low->a + s2*high->a;
  b = s1*low->b + s2*high->b;
  c = s1*low->c + s2*high->c;
  d = s1*low->d + s2*high->d;
  n = high->n;
  a -= R.a; b -= R.b; c -= R.c; d -= R.d; n -= R.n;
}
Elt & Elt::operator*= (const Elt & R) {
  int ax = a, bx = b, cx = c, dx = d;
  a = ax*R.a - bx*R.d - cx*R.c - dx*R.b;
  b = ax*R.b + bx*R.a - cx*R.d - dx*R.c;
  c = ax*R.c + bx*R.b + cx*R.a - dx*R.d;
  d = ax*R.d + bx*R.c + cx*R.b + dx*R.a;
  n += R.n;
}
const Elt Elt::operator+  (const Elt & R) const {
  Elt ret = *this;
  ret += R;
  return ret;
}
const Elt Elt::operator-  (const Elt & R) const {
  Elt ret = *this;
  ret -= R;
  return ret;
}
const Elt Elt::operator*  (const Elt & R) const {
  Elt ret = *this;
  ret *= R;
  return ret;
}

complex<double> Elt::to_complex() const {
  double e1, e2, e3;
  e1 = PI/4;
  e2 = PI/2;
  e3 = 3*PI/4;
  complex<double> ret(a + b*cos(e1) + c*cos(e2) + d*cos(e3),
                      b*sin(e1) + c*sin(e2) + d*sin(e3));
  return ret/pow(2, n);
}

void Elt::print() const {
  cout << "(" << a << " + " << b << "*w + " 
       << c << "*w^2 + " << d << "*w^3)/2^" << n;
}
