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

    Elt & operator=  (const Elt & R);
    Elt & operator+= (const Elt & R);
    Elt & operator-= (const Elt & R);
    Elt & operator*= (const Elt & R);
    const Elt operator+  (const Elt & R) const;
    const Elt operator-  (const Elt & R) const;
    const Elt operator*  (const Elt & R) const;

    complex<double> to_complex() const;
    void print() const;
};

