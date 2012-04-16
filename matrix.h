#ifndef MATRIX
#define MATRIX

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

    int rows() const;
    int cols() const;
    void resize(int m, int n);
   
    Rmatrix & operator=  (const Rmatrix & M);
    Rmatrix & operator+= (const Rmatrix & M);
    Rmatrix & operator-= (const Rmatrix & M);
    Rmatrix & operator*= (const Elt & R);
    Rmatrix & operator*= (const Rmatrix & M);
    Rmatrix & left_multiply (const Rmatrix & M);
    const Rmatrix operator+ (const Rmatrix & M) const;
    const Rmatrix operator- (const Rmatrix & M) const;
    const Rmatrix operator* (const Elt & R) const;
    const Rmatrix operator* (const Rmatrix & M) const;
    const bool operator== (const Rmatrix & M) const;
    Elt & operator() (int i, int j);

    const bool phase_eq(const Rmatrix & M) const;
    Unitary to_Unitary() const;
    void to_Unitary(Unitary & U) const;
    void adj(Rmatrix & M) const;
    void print() const;
    void submatrix(int m, int n, int numrow, int numcol, Rmatrix & M) const;
    bool is_nonlinear_reversible() const;
};

Rmatrix zero(int m, int n);
Rmatrix eye(int m, int n);

void matrix_test();

#endif
