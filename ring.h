#ifndef RING
#define RING

#include "configs.h"
#include <complex>
#include <iostream>
#include <stdlib.h>
#include <limits.h>

static const unsigned int reduce_bits[] = { 0xAAAAAAAA, 0xCCCCCCCC, 0xF0F0F0F0, 0xFF00FF00, 0xFFFF0000 };
static const int reduce_shift = sizeof(int) * CHAR_BIT - 1;
//////////////////

class Elt {
  private:
    int a, b, c, d, n;
  public:

    Elt::Elt() {a = b = c = d = n = 0;}
    Elt::Elt(int aa, int bb, int cc, int dd, int nn) { 
      a = aa; b = bb; c = cc; d = dd; n = nn;
    }
    Elt::Elt(const Elt & R) {
      *this = R;
    }
    static Elt randelt() { 
      return Elt(rand(), rand(), rand(), rand(), rand());
    }

    inline void reduce() {
      if (this->is_zero()) {
        n = 0;
        return;
      }
      int am = a >> reduce_shift;
      int bm = b >> reduce_shift;
      int cm = c >> reduce_shift;
      int dm = d >> reduce_shift;
      unsigned int tmp = ((a ^ am) - am) | ((b ^ bm) - bm) | ((c ^ cm) - cm) | ((d ^ dm) - dm);

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
    inline Elt & operator=  (const Elt & R) { a = R.a; b = R.b; c = R.c; d = R.d; n = R.n; }
    Elt & operator+= (const Elt & R);
    Elt & operator-= (const Elt & R);
    Elt & operator*= (const Elt & R);

    inline const Elt operator+  (const Elt & R) const {
      Elt ret = *this;
      ret += R;
      return ret;
    }
    inline const Elt operator-  (const Elt & R) const {
      Elt ret = *this;
      ret -= R;
      return ret;
    }
    inline const Elt operator*  (const Elt & R) const {
      Elt ret = *this;
      ret *= R;
      return ret;
    }
    inline const bool operator== (const Elt & R) const {
      return (a == R.a && b == R.b && c == R.c && d == R.d && n == R.n);
    }
    inline const bool operator!= (const Elt & R) const {
      return !(*this == R);
    }
    inline const bool operator<  (const Elt & R) const {
      if (a < R.a) return true;
      else if (a > R.a) return false;
      else if (b < R.b) return true;
      else if (b > R.b) return false;
      else if (c < R.c) return true;
      else if (c > R.c) return false;
      else if (d < R.d) return true;
      else if (d > R.d) return false;
      else if (n < R.n) return true;
      return false;
    }

    inline complex<double> to_complex() const {
      double rt = 1/sqrt(2);
      complex<double> ret(a + b*rt - d*rt, b*rt + c + d*rt);
      return ret/(double)(1 << n);
    }
    inline double abs() const {
      double rt = sqrt(2);
      double ret = a*a+b*b+c*c+d*d;
      ret += a*b*rt - a*d*rt + c*b*rt + c*d*rt;
      return sqrt(ret)/(1 << n);
    }
    inline Elt conj() { return Elt(a, -d, -c, -b, n); }
    inline bool is_zero() const {
      return (a == 0 && b == 0 && c == 0 && d == 0);
    }
    void print() const;

    friend const unsigned int hashfn(const pair<Elt, Elt> p);
    friend const bool eqfn(const pair<Elt, Elt> p, const pair<Elt, Elt> q);
};

struct elt_hasher {
  unsigned int operator()(const pair<Elt, Elt> p) const;
};
struct elt_eq {
  bool operator()(const pair<Elt, Elt> p, const pair<Elt, Elt> q) const;
};

void init_ring();

void test_ring();

#endif
