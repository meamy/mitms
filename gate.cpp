#define LA_COMPLEX_SUPPORT

#include <iostream>
#include <string>
#include <assert.h>
#include <map>

#include <gmc.h>
#include <blas3pp.h>
#include <blas2pp.h>
#include <blaspp.h>
#include <laslv.h>
#include <lavc.h>

/* Number of basis gates */
#define basis_size 9

/* Clifford group on 1 qubit */
#define I    0
#define H    1
#define X    2
#define Y    3
#define Z    4
#define S    5
#define Sd   6
#define T    7
#define Td   8

/* Basis state projections on one qubit */
#define PROJ(x, y)     (basis_size + 2*x + y)
#define IS_PROJ(x)  (basis_size <= x && x <= 0xfe)
#define GET_PROJ(x) (x - basis_size)

/* Controls. Highest order bit defines a control, the rest of the
   byte specifies the target */
#define C(x)          (x | 0x80)
#define IS_C(x)       (x & 0x80)
#define GET_TARGET(x) (x & 0x7F)

using namespace std;

int num_qubits = 0;
int dim = 0;

/* -------------- Type Definitions */
typedef LaGenMatComplex Unitary;

typedef struct {
  char qubits[2];
} Gate;

typedef struct Circuit {
	Gate G;
	struct Circuit * next;
} Circuit;
/*---------------------------------*/

const string gate_names[] = {
  "I",
  "H",
  "X",  
  "Y",
  "Z",  
  "S",
  "S*",
  "T",
  "T*",
};

const char adj[] = {
  I,
  H,
  X,
  Y,
  Z,
  Sd,
  S,
  Td,
  T,
};


Unitary * basis;

int fac(int n) {
  int ret = 1, i;
  for(i = n; i>1; i--) ret *= i;
  return ret;
}

void init(int n) {
  num_qubits = n;
  dim = pow(2, n);
  basis = new Unitary[basis_size + dim];

  int i;
  double rt = 1.0/sqrt(2.0);

	for (i = 1; i < basis_size + dim; i++) {
	  basis[i] = Unitary::zeros(2);
	}

  basis[I] = Unitary::eye(2);

	basis[H](0, 0) = LaComplex(rt, 0);
	basis[H](1, 0) = LaComplex(rt, 0);
	basis[H](0, 1) = LaComplex(rt, 0);
	basis[H](1, 1) = LaComplex(-rt, 0);

	basis[X](0, 1) = LaComplex(1, 0);
	basis[X](1, 0) = LaComplex(1, 0);

	basis[Y](0, 1) = LaComplex(0, -1);
	basis[Y](1, 0) = LaComplex(0, 1);

	basis[Z](0, 0) = LaComplex(1, 0);
	basis[Z](1, 1) = LaComplex(-1, 0);

	basis[S](0, 0) = LaComplex(1, 0);
	basis[S](1, 1) = LaComplex(0, 1);

  Blas_Mat_Mat_Mult(basis[I],   basis[S],   basis[Sd], false, true, 1, 0);

	basis[T](0, 0) = LaComplex(1, 0);
	basis[T](1, 1) = LaComplex(rt, rt);

  Blas_Mat_Mat_Mult(basis[I],   basis[T],   basis[Td], false, true, 1, 0);

  for(i = 0; i < dim; i++) {
    basis[PROJ((i / 2), (i % 2))](i%2, i%2) = LaComplex(1, 0);
  }
}

void print_gate(const Gate & G) {
  int i;
  char tmp;

  cout << "(";

  for (i = 0; i < num_qubits; i++) {
    tmp = G.qubits[i];

    if (IS_C(tmp)) {
      cout << "C(" << (int)GET_TARGET(tmp) << ")";
    } else {
      assert(tmp <= basis_size);
      cout << gate_names[tmp];
    }

    if (i < num_qubits - 1) { 
      cout << ".";
    } else {
      cout << ")";
    }
  }
}

void print_circuit(const Circuit * C) {
  const Circuit * tmp = C;
  while (tmp) {
    print_gate(tmp->G);
    cout << "\n";
    tmp = tmp->next;
  }
}

/* Compute the adjoint */
void gate_adj(const Gate & G, Gate & M) {
  int i;

  for (i = 0; i < num_qubits; i++) {
    if (IS_C(G.qubits[i])) {
      M.qubits[i] = G.qubits[i];
    } else {
      M.qubits[i] = adj[G.qubits[i]];
    }
  }
}

Circuit * circuit_adj(const Circuit * C, Circuit * last) {
  if (C == NULL) return last;

  Circuit * tmp = new Circuit;
  tmp->next = last;
  gate_adj(C->G, tmp->G);

  return circuit_adj(C->next, tmp);
}

/* Return the tensor product of the 1-qubit matrices */
Unitary tensor(const Gate & G) {
  const char * gates = G.qubits;
  int i, j, k, x, y;
  LaComplex tmp;

  assert(gates != NULL);
  for (i = 0; i < num_qubits; i++) {
    assert(gates[i] < basis_size + dim);
  }

  Unitary U(dim, dim);

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      tmp = LaComplex(1, 0);
      x = i; y = j;
      for (k = num_qubits-1; k >= 0; k--) {
        tmp *= (LaComplex)(basis[gates[k]](x % 2, y % 2));
        x /= 2; y /= 2;
      }
      U(i, j) = tmp;
    }
  }

  return U;
}

/* Compute a unitary for the given gate */
Unitary Unitary_of_Gate(const Gate & G);
Unitary Unitary_of_Gate(const Gate & G) {
  int i;
  const char *gates = G.qubits;

  for (i = 0; i < num_qubits; i++) {
    if (IS_C(gates[i])) {
      Gate A, B;
      char *a = A.qubits, *b = B.qubits;
      char tmp = GET_TARGET(gates[i]);

      for (int j = 0; j < num_qubits; j++) {
        if (j == i) {
          a[j] = PROJ(i, 0);
          b[j] = PROJ(i, 1);
        } else if (j == tmp) {
          a[j] = I;
          b[j] = gates[j];
        } else {
          a[j] = gates[j];
          b[j] = gates[j];
        }
      }

      return Unitary_of_Gate(A) + Unitary_of_Gate(B);
    }
  }
  return tensor(G);
}

Unitary Unitary_of_Circuit(const Circuit * C);
Unitary Unitary_of_Circuit(const Circuit * C) {
  assert(C != NULL);

  Unitary U = Unitary_of_Gate(C->G);
  Unitary V(dim, dim);
  if (C->next != NULL) {
    Blas_Mat_Mat_Mult(U, Unitary_of_Circuit(C->next), V, false, false, 1, 0);
    return V;
  } else {
    return U;
  }
}

double spec_norm(const Unitary & U) {
  Unitary G(U.size(0), U.size(1));
  Unitary F(U.size(0), U.size(1));
  LaVectorComplex e(G.size(1));
  double max = 0, temp;

  Blas_Mat_Mat_Mult(U, U, G, true, false, 1, 0);
  LaEigSolve(G, e, F);
 
  for (int i = 0; i < e.size(); i++) {
    temp = abs((complex<double>)((LaComplex)e(i)));
    temp = sqrt(temp);
    if (temp >= max) max = temp;
  }

  return max;
}

double dist(const Unitary & U, const Unitary & V) {
  #ifdef PHASE
    return spec_norm(U - V);
  #else
    double x, y, z;
    Unitary A(U.size(0), U.size(0));
    Unitary B(V.size(0), V.size(0));

    Blas_Mat_Mat_Mult(basis[X], U, A, false, true, 1, 0);
    Blas_Mat_Mat_Mult(U, A, A, false, false, 1, 0);
    Blas_Mat_Mat_Mult(basis[X], V, B, false, true, 1, 0);
    Blas_Mat_Mat_Mult(V, B, B, false, false, 1, 0);
    x = spec_norm(A - B);

    Blas_Mat_Mat_Mult(basis[Y], U, A, false, true, 1, 0);
    Blas_Mat_Mat_Mult(U, A, A, false, false, 1, 0);
    Blas_Mat_Mat_Mult(basis[Y], V, B, false, true, 1, 0);
    Blas_Mat_Mat_Mult(V, B, B, false, false, 1, 0);
    y = spec_norm(A - B);

    Blas_Mat_Mat_Mult(basis[Z], U, A, false, true, 1, 0);
    Blas_Mat_Mat_Mult(U, A, A, false, false, 1, 0);
    Blas_Mat_Mat_Mult(basis[Z], V, B, false, true, 1, 0);
    Blas_Mat_Mat_Mult(V, B, B, false, false, 1, 0);
    z = spec_norm(A - B);

    return sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
  #endif
}

int Hash_Unitary(const Unitary &U) {
  return (int)(10000000*dist(U, LaGenMatComplex::eye(dim)));
}

int main() {
  init(2);

  Gate G;
  G.qubits[0] = C(1);
  G.qubits[1] = H;
  print_gate(G);
  cout << "\n" << Unitary_of_Gate(G);

  return 0;
}
