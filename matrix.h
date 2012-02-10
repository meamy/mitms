#define LA_COMPLEX_SUPPORT

#include "ring.h"
#include <gmc.h>

typedef LaGenMatComplex Unitary;

class Rmatrix {
  private:
    int m, n;
    Elt ** mat;
  public:
	  Rmatrix();
		Rmatrix(int a, int b);
		Rmatrix(const Rmatrix & M);
	  ~Rmatrix();
   
    Rmatrix & operator=  (const Rmatrix & M);
    Rmatrix & operator+= (const Rmatrix & M);
    Rmatrix & operator-= (const Rmatrix & M);
    Rmatrix & operator*= (const Elt & R);
    Rmatrix & operator*= (const Rmatrix & M);
    const Rmatrix operator+ (const Rmatrix & M) const;
    const Rmatrix operator- (const Rmatrix & M) const;
    const Rmatrix operator* (const Elt & R) const;
    const Rmatrix operator* (const Rmatrix & M) const;

    Unitary to_Unitary() const;
    void print() const;
};

