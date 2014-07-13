
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

#ifndef MATRIX
#define MATRIX

#define LA_COMPLEX_SUPPORT

#include "ring.h"
#include <gmc.h>

typedef LaGenMatComplex Unitary;

class Rmatrix {
  private:
    int m, n;
    Elt * mat;
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

		Elt trace() const;
    const bool phase_eq(const Rmatrix & M) const;
    const bool equiv(const Rmatrix & M) const;

    Unitary to_Unitary() const;
    LaGenMatDouble to_Unitary_abs() const;
    Unitary to_Unitary_canon() const;
    void to_Unitary(Unitary & U) const;
    void to_Unitary_abs(LaGenMatDouble & U) const;
    void to_Unitary_canon(Unitary & U) const;

    void adj(Rmatrix & M) const;
    void permute(Rmatrix & M, int x) const;
    void permute_adj(Rmatrix & M, int x) const;
    void submatrix(int m, int n, int numrow, int numcol, Rmatrix & M) const;

    void print() const;
    void canon_phase();
    bool is_nonlinear_reversible() const;
};

Rmatrix zero(int m, int n);
Rmatrix eye(int m, int n);
Rmatrix col_permutation(int m, int n, int perm);

int fac(int n);
const char * from_lexi(int n);
int to_lexi(int i);

void init_rmatrix();
void test_rmatrix();

/* So that we can have functions that work for Rmatrixes and Unitaries */
void adj(const Rmatrix & M, Rmatrix & N);
void permute(const Rmatrix & M, Rmatrix & N, int x);
void permute_adj(const Rmatrix & M, Rmatrix & N, int x);
void submatrix(const Rmatrix & M, int m, int n, int numrow, int numcol, Rmatrix & N);

void adj(const Unitary & M, Unitary & N);
void permute(const Unitary & M, Unitary & N, int x);
void permute_adj(const Unitary & M, Unitary & N, int x);
//void submatrix(const Unitary & M, int m, int n, int numrow, int numcol, Unitary & N);
Unitary operator*(const Unitary & U, const Unitary & V);

#endif
