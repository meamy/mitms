
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

#ifndef RING
#define RING

#include "configs.h"
#include <complex>
#include <iostream>
#include <stdlib.h>
#include <limits.h>

static const double rt_two = 1/sqrt(2);
//////////////////

class Elt {
  private:
    int a, b, c, d, n;
  public:

    inline Elt() {a = b = c = d = n = 0;}
    inline Elt(int aa, int bb, int cc, int dd, int nn) { 
      a = aa; b = bb; c = cc; d = dd; n = nn;
    }
    inline Elt(const Elt & R) {
      *this = R;
    }

    static Elt randelt() { 
      return Elt(rand(), rand(), rand(), rand(), rand());
    }

    /* Arithmetic operations -- inlining gave worse performance somehow */
    void reduce();
    inline Elt & operator=  (const Elt & R) {
      a = R.a; b = R.b; c = R.c; d = R.d; n = R.n;
    }
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
		/*
    inline Elt            norm() const {
			int tmp = a*(b - d) + c*(b + d);
      return Elt(a*a + b*b + c*c + d*d, tmp, 0, -tmp, 2*n);
    }
		*/
    inline double         norm() const {
			if (n < 0) {
				return (a*a + b*b + c*c + d*d + (a*(b-d) + c*(b+d))*rt_two) * (double)(1 << 2*(-n));
			}else {
				return (a*a + b*b + c*c + d*d + (a*(b-d) + c*(b+d))*rt_two) / (double)(1 << 2*n);
			}
    }
    inline double          abs() const {
      return sqrt(a*a + b*b + c*c + d*d + (a*(b-d) + c*(b+d))*rt_two) / (double)(1 << n);
    }
    inline Elt            conj() const { 
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
