#ifndef RING
#define RING

#include "configs.h"
#include <complex>
#include <iostream>
#include <stdlib.h>
#include <limits.h>

static const unsigned int reduce_bits[] = { 0xAAAAAAAA, 0xCCCCCCCC, 0xF0F0F0F0, 0xFF00FF00, 0xFFFF0000 };
static const int reduce_shift = sizeof(int) * CHAR_BIT - 1;
static const double rt_two = 1/sqrt(2);
//////////////////

class Elt {
  private:
    int a, b, c, d, n;
  public:

    Elt(); 
    Elt(int aa, int bb, int cc, int dd, int nn);
    Elt(const Elt & R);
    static Elt randelt() { 
      return Elt(rand(), rand(), rand(), rand(), rand());
    }

    /* Arithmetic operations */
    void reduce();
    Elt & operator=  (const Elt & R);
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

    /* Comparison functions */
    inline const bool is_zero() const { 
      return (a == 0 && b == 0 && c == 0 && d == 0);
    }
    inline const bool operator== (const Elt & R) const {
      return (a == R.a && b == R.b && c == R.c && d == R.d && n == R.n);
    }
    inline const bool operator!= (const Elt & R) const {
      return (a != R.a || b != R.b || c != R.c || d != R.d || n != R.n);
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

    /* Transformations */
    inline complex<double> to_complex() const {
      return complex<double>(a + (b-d)*rt_two, (b+d)*rt_two + c) / (double)(1 << n);
    }
    inline double abs() const {
      return sqrt(a*a + b*b + c*c + d*d + (a*(b-d) + c*(b+d))*rt_two) / (double)(1 << n);
    }
    inline Elt conj() { 
      return Elt(a, -d, -c, -b, n);
    }

    void print() const;

    /* Hashing related */
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
