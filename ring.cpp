#include "ring.h"
#include <assert.h>
#include <unordered_map>
#include <limits.h>

#define PI M_PI

//////////////////
typedef pair<Elt, Elt> eltpair;
typedef unordered_map<eltpair, Elt, elt_hasher, elt_eq> hashmap;
hashmap mult_table;
//////////////////
static const unsigned int bits[] = { 0xAAAAAAAA, 0xCCCCCCCC, 0xF0F0F0F0, 0xFF00FF00, 0xFFFF0000 };
int shift = sizeof(int) * CHAR_BIT - 1;
//////////////////

Elt::Elt() {a = b = c = d = n = 0;}
Elt::Elt(int aa, int bb, int cc, int dd, int nn) { 
  a = aa; b = bb; c = cc; d = dd; n = nn;
}
Elt::Elt(const Elt & R) {
  *this = R;
}

void Elt::reduce() {
  if (this->is_zero()) {
    n = 0;
    return;
  }
  int am = a >> shift;
  int bm = b >> shift;
  int cm = c >> shift;
  int dm = d >> shift;
  int tmp = ((a ^ am) - am) | ((b ^ bm) - bm) | ((c ^ cm) - cm) | ((d ^ dm) - dm);
  int x = tmp & (~tmp + 1);
  if (x > 1) {
    tmp = (x & bits[0]) != 0;
    for (int i = 4; i > 0; i--) {
      tmp |= ((x & bits[i]) != 0) << i;
    }

    a /= x;
    b /= x;
    c /= x;
    d /= x;
    n -= tmp;
  }
/*
  int i, x;
  for (i = n; i > 0; i--) {
    x = pow(2, i);

    if (a % x == 0 && b % x == 0 && c % x == 0 && d % x == 0) {
      a /= x;
      b /= x;
      c /= x;
      d /= x;
      n -= i;
      return;
    } 
  }
  */
}

Elt & Elt::operator=  (const Elt & R) {
  a = R.a; b = R.b; c = R.c; d = R.d; n = R.n;
}
Elt & Elt::operator+= (const Elt & R) {
  /*
  const Elt * low, * high;
  int diff;
  if (n < R.n) {
    low = this;
    high = & R;
    diff = 1 << (R.n - n);
  } else {
    low = & R;
    high = this;
    diff = 1 << (n - R.n);
  }

  a = (diff*low->a) + high->a;
  b = (diff*low->b) + high->b;
  c = (diff*low->c) + high->c;
  d = (diff*low->d) + high->d;
  n = high->n;
  this->reduce();
  */
  int x = 1 << R.n;
  int y = 1 << n;
  n += R.n;
  a = x*a + y*R.a;
  b = x*b + y*R.b;
  c = x*c + y*R.c;
  d = x*d + y*R.d;
  this->reduce();
}
Elt & Elt::operator-= (const Elt & R) {
  /*
  const Elt * low, * high;
  int s1, s2;
  if (n < R.n) {
    low = this;
    high = & R;
    s1 = 1 << (R.n - n);
    s2 = -1;
  } else {
    low = & R;
    high = this;
    s1 = -1 << (n - R.n);
    s2 = 1;
  }

  a = s1*low->a + s2*high->a;
  b = s1*low->b + s2*high->b;
  c = s1*low->c + s2*high->c;
  d = s1*low->d + s2*high->d;
  n = high->n;
  a -= R.a; b -= R.b; c -= R.c; d -= R.d; n -= R.n;
  this->reduce();
  */
  int x = 1 << R.n;
  int y = 1 << n;
  n += R.n;
  a = x*a - y*R.a;
  b = x*b - y*R.b;
  c = x*c - y*R.c;
  d = x*d - y*R.d;
  this->reduce();
}
Elt & Elt::operator*= (const Elt & R) {
#if HASH
  if (mult_table.bucket_count() < 5) {
    mult_table.reserve(100000);
  }
  hashmap::iterator it = mult_table.find(eltpair(*this, R));
  if (it != mult_table.end()) {
    *this = it->second;
  } else {
#endif
    int ax = a, bx = b, cx = c, dx = d;
    Elt tmp = *this;
    a = ax*R.a - bx*R.d - cx*R.c - dx*R.b;
    b = ax*R.b + bx*R.a - cx*R.d - dx*R.c;
    c = ax*R.c + bx*R.b + cx*R.a - dx*R.d;
    d = ax*R.d + bx*R.c + cx*R.b + dx*R.a;
    n += R.n;
    this->reduce();
#if HASH
    mult_table.insert(pair<eltpair, Elt>(eltpair(tmp, R), *this));
  }
#endif
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

const bool Elt::operator<  (const Elt & R) const {
  if (a < R.a) return true;
  if (a > R.a) return false;
  if (b < R.b) return true;
  if (b > R.b) return false;
  if (c < R.c) return true;
  if (c > R.c) return false;
  if (d < R.d) return true;
  if (d > R.d) return false;
  if (n < R.n) return true;
  return false;
}

complex<double> Elt::to_complex() const {
  double rt = 1/sqrt(2);
  complex<double> ret(a + b*rt - d*rt, b*rt + c + d*rt);
  return ret/(double)(1 << n);
}

double Elt::abs() const {
  double rt = sqrt(2);
  double ret = pow(a,2)+pow(b,2)+pow(c,2)+pow(d,2);
  ret += a*b*rt - a*d*rt + c*b*rt + c*d*rt;
  return sqrt(ret)/(1 << n);
}

Elt Elt::conj() {
  return Elt(a, -d, -c, -b, n);
}

bool Elt::is_zero() const {
  return (a == 0 && b == 0 && c == 0 && d == 0);
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

const unsigned int hashfn(const pair<Elt, Elt> p) {
  return 
    abs(p.first.a + p.second.a) +
    abs(p.first.b + p.second.b)*10 +
    abs(p.first.c + p.second.c)*100 +
    abs(p.first.d + p.second.d)*1000 +
    abs(p.first.n + p.second.n)*10000;
}

const bool eqfn(const pair<Elt, Elt> p, const pair<Elt, Elt> q) {
  return ((p.first == q.first) && (p.second == q.second)) ||
         ((p.first == q.second) && (p.second == q.first));
}

unsigned int elt_hasher::operator()(const pair<Elt, Elt> p) const {
  return hashfn(p);
}

bool elt_eq::operator()(const pair<Elt, Elt> p, const pair<Elt, Elt> q) const {
  return eqfn(p, q);
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

