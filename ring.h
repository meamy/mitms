#ifndef RING
#define RING

#include "configs.h"
#include <complex>
#include <iostream>
#include <stdlib.h>

class Elt {
  private:
    int a, b, c, d, n;
  public:

    Elt();
    Elt(int a, int b, int c, int d, int n);
    Elt(const Elt & R);
    static Elt randelt() { 
      return Elt(rand(), rand(), rand(), rand(), rand());
    }

    void reduce();
    Elt & operator=  (const Elt & R);
    Elt & operator+= (const Elt & R);
    Elt & operator-= (const Elt & R);
    Elt & operator*= (const Elt & R);
    const Elt operator+  (const Elt & R) const;
    const Elt operator-  (const Elt & R) const;
    const Elt operator*  (const Elt & R) const;
    const bool operator== (const Elt & R) const;
    const bool operator!= (const Elt & R) const;
    const bool operator<  (const Elt & R) const;

    complex<double> to_complex() const;
    double abs() const;
    Elt conj();
    bool is_zero() const;
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
