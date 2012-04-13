#include "ring.h"
#include <assert.h>

#define PI M_PI

Elt::Elt() {a = b = c = d = n = 0;}
Elt::Elt(int aa, int bb, int cc, int dd, int nn) { 
  a = aa; b = bb; c = cc; d = dd; n = nn;
}
Elt::Elt(const Elt & R) {
  *this = R;
}

void Elt::reduce() {
  int i, x;
  for (i = n; i > 0; i--) {
    x = pow(2, i);
    if (a % x == 0 && b % x == 0 && c % x == 0 && d % x == 0) {
      a /= x;
      b /= x;
      c /= x;
      d /= x;
      n -= i;
      break;
    }
  }
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
  this->reduce();
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
  this->reduce();
}
Elt & Elt::operator*= (const Elt & R) {
  int ax = a, bx = b, cx = c, dx = d;
  a = ax*R.a - bx*R.d - cx*R.c - dx*R.b;
  b = ax*R.b + bx*R.a - cx*R.d - dx*R.c;
  c = ax*R.c + bx*R.b + cx*R.a - dx*R.d;
  d = ax*R.d + bx*R.c + cx*R.b + dx*R.a;
  n += R.n;
  this->reduce();
}
const Elt Elt::operator+  (const Elt & R) const {
  Elt ret = *this;
  ret += R;
  ret.reduce();
  return ret;
}
const Elt Elt::operator-  (const Elt & R) const {
  Elt ret = *this;
  ret -= R;
  ret.reduce();
  return ret;
}
const Elt Elt::operator*  (const Elt & R) const {
  Elt ret = *this;
  ret *= R;
  ret.reduce();
  return ret;
}

const bool Elt::operator== (const Elt & R) const {
  return (a == R.a && b == R.b && c == R.c && d == R.d && n == R.n);
}

const bool Elt::operator!= (const Elt & R) const {
  return !(*this == R);
}

complex<double> Elt::to_complex() const {
  double rt = 1/sqrt(2);
  complex<double> ret(a + b*rt - d*rt, b*rt + c + d*rt);
  return ret/pow(2, n);
}

Elt Elt::conj() {
  return Elt(a, -d, -c, -b, n);
}

void Elt::print() const {
  if (a == 0 && b == 0 && c == 0 && d == 0) {
    cout << "0";
  } else {
    if (n != 0) cout << "(";
    if (a != 0) {
      cout << a;
      if (b != 0 || c != 0 || d != 0) cout << " + ";
    }
    if (b != 0) {
      cout << b << "w";
      if (c != 0 || d != 0) cout << " + ";
    }
    if (c != 0) {
      cout << c << "w^2";
      if (d != 0) cout << " + ";
    }
    if (d != 0) cout << d << "w^3";
    if (n != 0) cout << ")/2^" << n;
  }
}

void ring_test() {
  Elt a;
  Elt b(1, 1, 1, 1, 2);
  Elt c(b);

  assert(b == c);
  c *= Elt(0, 0, 1, 0, 4);
  assert(b != c);
  a = b.conj();
  c = a + b;
  assert(imag(c.to_complex()) == 0);
  c = a - b;
  assert(b == Elt(1, 1, 1, 1, 2));
}

