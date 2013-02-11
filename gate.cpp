#include "gate.h"
#include <iomanip>

bool identities[basis_size][basis_size];

typedef Rmatrix* Rptr;
const Rptr * gate_ht;
const Rmatrix * basis;
unsigned int max_hash = 0;

/* -------------- Gates */
/* Non-exported functions */
bool valid_gate(gate G) {
  int i, j, y;
  char x;
  bool flg;
  for (i = 0; i < num_qubits; i++) {
    x = G[i];
    if (x == X || x == Y || x == Z) {
      if (!config::paulis) {
        flg = false;
        for(j = 0; j < num_qubits; j++) {
          flg = flg || (IS_C(G[j]) && GET_TARGET(G[j]) == i);
        }
        if (flg == false) return false;
      }
    } else if (IS_C(x)) {
      y = GET_TARGET(x);
      if (y == -1 || G[y] != X) return false;
      for (j = i+1; j < num_qubits; j++) {
        if (IS_C(G[j]) && GET_TARGET(G[j]) == y) return false;
      }
    }
  }
  return true;
}


void increment(gate G, bool t) {
  int i, j, x;
  for (i = num_qubits-1; i >= 0; i--) {
      /* gate i is a single qubit gate, and the next gate is also single qubit */
    if (0 <= G[i]  && G[i] < (basis_size - 1 - 2*((int)(!t)))) {
      G[i] = G[i] + 1;
      return;
      /* the next gate is a CNOT -- make it an invalid control to be coupled with X's later */
    } else if (G[i] == basis_size - 1 - 2*((int)(!t))) {
      G[i] = C(0);
      return;
      /* the gate is a CNOT and we've gone through all combinations of CNOTs on qubits after i */
    } else if (IS_C(G[i])) {
      j = GET_TARGET(G[i]);
      if (j < num_qubits - 1) {
        G[i] = C(j+1);
        return;
      } else {
        G[i] = I;
      }
    } else {
      assert(false);
    }
  }
  return;
}

/* Return the tensor product of the 1-qubit matrices */
void tensor(const gate G, Rmatrix & U, bool adj) {
  int i, j, k, x, y;
  Elt tmp;

  for (i = 0; i < num_qubits; i++) {
    assert(G[i] < basis_size + dim);
  }

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      tmp = Elt(1, 0, 0, 0, 0);
      x = i; y = j;
      for (k = num_qubits-1; k >= 0; k--) {
        tmp *= basis[G[k]](x % 2, y % 2);
        x /= 2; y /= 2;
      }
      if (adj) {
        U(j, i) = tmp.conj();
        U(j, i).reduce();
      } else {
        U(i, j) = tmp;
        U(i, j).reduce();
      }
    }
  }
}

inline void tensor(const gate G, Rmatrix & U) {
  tensor(G, U, false);
}

/* Exported functions */
void next_gate(gate G) {
  increment(G, true);
  while(!valid_gate(G)) increment(G, true);
}

void next_clifford(gate G) {
  increment(G, false);
  while(!valid_gate(G)) increment(G, false);
}

bool gate_eq(const gate A, const gate B) {
  int i;
  for (i = 0; i < num_qubits; i++) {
    if (A[i] != B[i]) return false;
  }
  return true;
}

bool is_eye(const gate G) {
  int i;
  for (i = 0; i < num_qubits; i++) {
    if (G[i] != I) return false;
  }
  return true;
}

void gate_transform(const gate A, gate B, const char * perm, bool adj) {
  char gt;

  for (int i = 0; i < num_qubits; i++) {
    gt = A[i];
    if (IS_C(gt)) {
      B[perm[i]] = C(perm[GET_TARGET(gt)]);
    } else {
      B[perm[i]] = adj ? adjoint[gt] : gt;
    }
  }
}

void gate_transform(const gate A, gate B, int i, bool adj) {
  const char * perm = from_lexi(i);
  gate_transform(A, B, perm, adj);
}

void expand_cnots(const gate G, Rmatrix & U, bool adj);
void expand_cnots(const gate G, Rmatrix & U, bool adj) {
  int i;
  for (i = 0; i < num_qubits; i++) {
    if (IS_C(G[i])) {
      gate A, B;
      A = new char[num_qubits];
      B = new char[num_qubits];
      Rmatrix V(dim, dim);
      char tmp = GET_TARGET(G[i]);

      for (int j = 0; j < num_qubits; j++) {
        if (j == i) {
          A[j] = PROJ(i, 0);
          B[j] = PROJ(i, 1);
        } else if (j == tmp) {
          A[j] = I;
          B[j] = G[j];
        } else {
          A[j] = G[j];
          B[j] = G[j];
        }
      }

      expand_cnots(A, V, adj);
      expand_cnots(B, U, adj);
      U += V;
      delete [] A;
      delete [] B;

      return;
    }
  }
  tensor(G, U, adj);
}

/* Compute a unitary for the given gate */
void gate_to_Rmatrix(const gate G, Rmatrix & U, bool adj) {
  int i;
  if (!config::tensors) {
    if (adj) {
      gate A = new char[num_qubits];
      gate_transform(G, A, 0, true);
      i = gate_hasher(A);
      delete [] A;
    } else {
      i = gate_hasher(G);
    }
    if (i <= max_hash && gate_ht[i] != NULL) {
      U = *(gate_ht[i]);
      return;
    }
  }
  expand_cnots(G, U, adj);
}

void gate_to_Unitary(const gate G, Unitary & U, bool adj) {
  Rmatrix tmp(dim, dim);
  gate_to_Rmatrix(G, tmp, adj);
  tmp.to_Unitary(U);
}


void print_gate(gate G) {
  int i;
  char tmp;

  for (i = 0; i < num_qubits; i++) {
    tmp = G[i];

    if (IS_C(tmp)) {
      cout << "C(" << (int)GET_TARGET(tmp) + 1 << ")";
    } else {
      assert(tmp <= basis_size);
      cout << " " << setw(3) << left << gate_names[tmp];
    }

    cout << "\n";
  }
}

void output_gate(const gate G, ofstream & out) {
  out.write(G, num_qubits);
}
void input_gate(gate G, ifstream & in) {
  in.read(G, num_qubits);
}

/* Determine if there is a nontrivial pair of gates that multiply to the identity */
bool nontrivial_id(const gate A, const gate B) {
  int i, j, mask = ~((int)0), tmp = 0, x;
  bool ret = false;
  for (i = 0; i < num_qubits && !ret; i++) {
    if ((1 << i) & mask) {
      if (IS_C(B[i])) {
        x = ~(1 << GET_TARGET(B[i]));
        mask &= x;
        tmp &= x; 
        ret = ret || (A[i] == B[i]);
      } else if (IS_C(A[i])) {
        x = ~(1 << GET_TARGET(A[i]));
        mask &= x;
        tmp &= x; 
        ret = ret || (A[i] == B[i]);
      } else if (A[i] == X && B[i] == X) {
        tmp |= 1 << i;
      } else {
        ret = ret || identities[A[i]][B[i]];
      }
    }
  }
  return ret || (bool)tmp;
}

unsigned int gate_hasher(const gate G) {
  unsigned int ret = 0;
  unsigned int tmp = 0;
  int i;
  for (i = 0; i < num_qubits; i++) {
    ret *= basis_size + num_qubits;
    if (IS_C(G[i])) {
      ret += basis_size + GET_TARGET(G[i]);
    } else if (IS_PROJ(G[i])) {
      return 0;
    } else {
      ret += G[i];
    }
  }
  return ret;
}

void init_identities() {
  int i, j;
  for (i = 0; i < basis_size; i++) {
    for (j = 0; j < basis_size; j++) {
      if (i == T && (j == Sd)) {
        identities[i][j] = true;
      } else if (i == Td && (j == S)) {
        identities[i][j] = true;
      } else if (i == Sd && j == T) {
        identities[i][j] = true;
      } else if (i == S && j == Td) {
        identities[i][j] = true;
      } else if (i == adjoint[j]) {
        identities[i][j] = true;
      } else {
        identities[i][j] = false;
      }
    }
  }
}

void init_gate() {
  /*-------------------------- Initializing matrices */
  int i;
  Rmatrix * basis_tmp = new Rmatrix[basis_size + dim];
	for (i = 1; i < basis_size + dim; i++) {
	  basis_tmp[i] = zero(2, 2);
	}

  basis_tmp[I] = eye(2, 2);

	basis_tmp[H](0, 0) = Elt(0, 1, 0, -1, 1);
	basis_tmp[H](1, 0) = Elt(0, 1, 0, -1, 1);
	basis_tmp[H](0, 1) = Elt(0, 1, 0, -1, 1);
	basis_tmp[H](1, 1) = Elt(0, -1, 0, 1, 1);

	basis_tmp[X](0, 1) = Elt(1, 0, 0, 0, 0);
	basis_tmp[X](1, 0) = Elt(1, 0, 0, 0, 0);

	basis_tmp[Y](0, 1) = Elt(0, 0, -1, 0, 0);
	basis_tmp[Y](1, 0) = Elt(0, 0, 1, 0, 0);

	basis_tmp[Z](0, 0) = Elt(1, 0, 0, 0, 0);
	basis_tmp[Z](1, 1) = Elt(-1, 0, 0, 0, 0);

	basis_tmp[S](0, 0) = Elt(1, 0, 0, 0, 0);
	basis_tmp[S](1, 1) = Elt(0, 0, 1, 0, 0);

  basis_tmp[S].adj(basis_tmp[Sd]);

	basis_tmp[T](0, 0) = Elt(1, 0, 0, 0, 0);
	basis_tmp[T](1, 1) = Elt(0, 1, 0, 0, 0);

  basis_tmp[T].adj(basis_tmp[Td]);

  for (i = 0; i < dim; i++) {
    basis_tmp[PROJ((i / 2), (i % 2))](i%2, i%2) = Elt(1, 0, 0, 0, 0);
  }

  basis = basis_tmp;

  /*--------------------------- Initializing hash table */
  if (!config::tensors) {
    config::tensors = true;
    gate G = new char[num_qubits];
    Rmatrix R(dim, dim);
    max_hash = (int)pow((double)(basis_size + num_qubits), num_qubits);
    Rptr * tmp = new Rptr [max_hash];

    int j;
    for (j = 0; j < max_hash; j++) tmp[j] = NULL;
    for (j = 0; j < num_qubits; j++) G[j] = I;

    gate_to_Rmatrix(G, R, false);
    tmp[gate_hasher(G)] = new Rmatrix(dim, dim);
    *(tmp[gate_hasher(G)]) = R;

    next_gate(G);
    while(!is_eye(G)) {
      gate_to_Rmatrix(G, R, false);
			j = gate_hasher(G);
			if (tmp[j] != NULL) {
				cout << "ERROR: hash collision at " << j << "\n";
				print_gate(G);
				exit(1);
			}
      tmp[j] = new Rmatrix(dim, dim);
      *(tmp[j]) = R;
      next_gate(G);
    }

    gate_ht = tmp;
    config::tensors = false;
  }

  /*--------------------------- Initializing the identitites */
  init_identities();
}


void test_gate() {
  gate A = new char[num_qubits];
  gate B = new char[num_qubits];
  gate C = new char[num_qubits];

  A[0] = H;
  A[1] = I;
  A[2] = X;

  copy_gate(A, B);
  assert(gate_eq(B, A));
  B[0] = I;
  assert(gate_neq(B, A));

  gate_adj(A, C);
  const char * tmp = from_lexi(2);
  gate_permute(A, B, tmp);

  Rmatrix R(dim, dim);
  Unitary U(dim, dim);
  gate_to_Rmatrix(B, R);
  gate_to_Unitary(C, U);

  delete [] A;
  delete [] B;
  delete [] C;
}

