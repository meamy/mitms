
/*--------------------------------------------------------------------
MITMS - Meet-in-the-middle quantum circuit synthesis
Copyright (C) 2013  Matthew Amy and The University of Waterloo,
Institute for Quantum Computing, Quantum Circuits Group

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Author: Matthew Amy
---------------------------------------------------------------------*/

#include "ring.h"
#include <assert.h>

static const unsigned int reduce_bits[] = { 0xAAAAAAAA, 0xCCCCCCCC, 0xF0F0F0F0, 0xFF00FF00, 0xFFFF0000 };
static const int reduce_shift = sizeof(int) * CHAR_BIT - 1;

typedef pair<Elt, Elt> eltpair;
typedef hash_table<eltpair, Elt, elt_hasher, elt_eq>::t hashmap;
hashmap * mult_table;

/* Arithmetic operations */
/* Reduction to lowest form is handled manually by programs using Elts rather than
   after every arithmetic operation. Users should be careful to reduce often enough
   to avoid integer overflow, but also minimize their usage of this function */
void Elt::reduce() {
  if (this->is_zero()) {
    n = 0;
    return;
  }
  int am = a >> reduce_shift;
  int bm = b >> reduce_shift;
  int cm = c >> reduce_shift;
  int dm = d >> reduce_shift;
  unsigned int tmp = ((a ^ am) - am) | 
    ((b ^ bm) - bm) | ((c ^ cm) - cm) | ((d ^ dm) - dm);

  int x = tmp & (~tmp + 1);
  if (x > 1) {
    tmp = (x & reduce_bits[0]) != 0;
    for (int i = 4; i > 0; i--) {
      tmp |= ((x & reduce_bits[i]) != 0) << i;
    }

    a = a >> tmp;
    b = b >> tmp;
    c = c >> tmp;
    d = d >> tmp;
    n -= tmp;
  }
}
Elt & Elt::operator+= (const Elt & R) {
  int x = 1 << R.n;
  int y = 1 << n;
  n += R.n;
  a = x*a + R.a*y;
  b = x*b + R.b*y;
  c = x*c + R.c*y;
  d = x*d + R.d*y;
}
Elt & Elt::operator-= (const Elt & R) {
  int x = 1 << R.n;
  int y = 1 << n;
  n += R.n;
  a = x*a - y*R.a;
  b = x*b - y*R.b;
  c = x*c - y*R.c;
  d = x*d - y*R.d;
}
Elt & Elt::operator*= (const Elt & R) {
  if (config::hash_ring) {
    hashmap::iterator it = mult_table->find(eltpair(*this, R));
    if (it != mult_table->end()) {
      *this = it->second;
    } else {
      int ax = a, bx = b, cx = c, dx = d;
      Elt tmp = *this;
      a = ax*R.a - bx*R.d - cx*R.c - dx*R.b;
      b = ax*R.b + bx*R.a - cx*R.d - dx*R.c;
      c = ax*R.c + bx*R.b + cx*R.a - dx*R.d;
      d = ax*R.d + bx*R.c + cx*R.b + dx*R.a;
      n += R.n;
      mult_table->insert(pair<eltpair, Elt>(eltpair(tmp, R), *this));
    }
  } else {
    int ax = a, bx = b, cx = c, dx = d;
    a = ax*R.a - bx*R.d - cx*R.c - dx*R.b;
    b = ax*R.b + bx*R.a - cx*R.d - dx*R.c;
    c = ax*R.c + bx*R.b + cx*R.a - dx*R.d;
    d = ax*R.d + bx*R.c + cx*R.b + dx*R.a;
    n += R.n;
  }
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

void init_ring() {
  if (config::hash_ring) {
    mult_table = new hashmap;
  }
}

void test_ring() {
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

