#include <complex>
#include <iostream>

using namespace std;

class Ring {
  private:
    int a, b, c, d, n;
  public:

    Ring();
    Ring(int a, int b, int c, int d, int n);
    Ring(const Ring & R);

    Ring & operator=  (const Ring & R);
    Ring & operator+= (const Ring & R);
    Ring & operator-= (const Ring & R);
    Ring & operator*= (const Ring & R);
    const Ring operator+  (const Ring & R) const;
    const Ring operator-  (const Ring & R) const;
    const Ring operator*  (const Ring & R) const;

    complex<double> to_complex() const;
    void print() const;
};

