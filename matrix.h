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

    static Rmatrix rand(int m, int n);

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
    const bool operator<  (const Rmatrix & M) const;
    const Elt & operator() (int i, int j) const;
    Elt & operator() (int i, int j);

    const bool phase_eq(const Rmatrix & M) const;
    Unitary to_Unitary() const;
    LaGenMatDouble to_Unitary_abs() const;
    Unitary to_Unitary_canon() const;
    void to_Unitary(Unitary & U) const;
    void to_Unitary_abs(LaGenMatDouble & U) const;
    void to_Unitary_canon(Unitary & U) const;
    void adj(Rmatrix & M) const;
    void permute(Rmatrix & M, int x) const;
    void permute_adj(Rmatrix & M, int x) const;
    void print() const;
    void submatrix(int m, int n, int numrow, int numcol, Rmatrix & M) const;
    void canon_phase();
    bool is_nonlinear_reversible() const;
};

Rmatrix zero(int m, int n);
Rmatrix eye(int m, int n);

int fac(int n);
char * from_lexi(int n);
int to_lexi(int i);
int inv_permutation(int i);
void init_permutations(int n);
void matrix_test();

#endif
