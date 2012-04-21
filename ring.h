#ifndef RING
#define RING

#include <complex>
#include <iostream>

using namespace std;

class Elt {
  private:
    int a, b, c, d, n;
  public:

    Elt();
    Elt(int a, int b, int c, int d, int n);
    Elt(const Elt & R);
    static Elt rand() { 
      return Elt(std::rand(), std::rand(), std::rand(), std::rand(), std::rand() % 5);
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
};

void ring_test();

#endif
