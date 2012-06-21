#include "gate.h"
#include <iomanip>

bool identities[basis_size][basis_size];

const Rmatrix * gate_ht;
const Rmatrix * basis;

/* -------------- Gates */
Gate & Gate::operator=(const Gate & G) {
  if (gates == NULL) {
    gates = new char[num_qubits];
  }
  for (int i = 0; i < num_qubits; i++) {
    gates[i] = G.gates[i];
  }
  return *this;
}

char & Gate::operator[](int i) const { 
  assert(0 <= i && i < num_qubits);
  return gates[i];
}

bool Gate::valid_gate() {
  int i, j, y;
  char x;
  bool flg;
  for (i = 0; i < num_qubits; i++) {
    x = gates[i];
    if (x == X || x == Y || x == Z) {
      if (!config::paulis) {
        flg = false;
        for(j = 0; j < num_qubits; j++) {
          flg = flg || (IS_C(gates[j]) && GET_TARGET(gates[j]) == i);
        }
        if (flg == false) return false;
      }
    } else if (IS_C(x)) {
      y = GET_TARGET(x);
      if (y == -1 || gates[y] != X) return false;
      for (j = i+1; j < num_qubits; j++) {
        if (IS_C(gates[j]) && GET_TARGET(gates[j]) == y) return false;
      }
    }
  }
  return true;
}


void Gate::increment(bool t) {
  int i, j, x;
  for (i = num_qubits-1; i >= 0; i--) {
      /* gate i is a single qubit gate, and the next gate is also single qubit */
    if (0 <= gates[i]  && gates[i] < (basis_size - 1 - 2*((int)(!t)))) {
      gates[i] = gates[i] + 1;
      return;
      /* the next gate is a CNOT -- make it an invalid control to be coupled with X's later */
    } else if (gates[i] == basis_size - 1 - 2*((int)(!t))) {
      gates[i] = C(0);
      return;
      /* the gate is a CNOT and we've gone through all combinations of CNOTs on qubits after i */
    } else if (IS_C(gates[i])) {
      j = GET_TARGET(gates[i]);
      if (j < num_qubits - 1) {
        gates[i] = C(j+1);
        return;
      }
      gates[i] = I;
    } else {
      assert(false);
    }
  }
  return;
}

Gate & Gate::operator++() {
  this->increment(true);
  while(!(this->valid_gate())) this->increment(true);
  return *this;
}

Gate & Gate::cliffpp() {
  this->increment(false);
  while(!(this->valid_gate())) this->increment(false);
  return *this;
}

const bool Gate::operator==(const Gate & G) const {
  int i;
  for (i = 0; i < num_qubits; i++) {
    if (gates[i] != G[i]) return false;
  }
  return true;
}

Gate::Gate() { gates = new char[num_qubits]; }
Gate::Gate(const Gate & G) { gates = new char[num_qubits]; *this = G; }
Gate::~Gate() { delete [] gates; }

const bool Gate::eye() const {
  int i;

  for (i = 0; i < num_qubits; i++) {
    if (gates[i] != I) return false;
  }

  return true;
}

void Gate::adj(Gate & G) const {
  int i;

  for (i = 0; i < num_qubits; i++) {
    if (IS_C(gates[i])) {
      G[i] = gates[i];
    } else {
      G[i] = adjoint[gates[i]];
    }
  }
}

/* Return the tensor product of the 1-qubit matrices */
void Gate::tensor(Rmatrix & U) const {
  int i, j, k, x, y;
  Elt tmp;

  for (i = 0; i < num_qubits; i++) {
    assert(gates[i] < basis_size + dim);
  }

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      tmp = Elt(1, 0, 0, 0, 0);
      x = i; y = j;
      for (k = num_qubits-1; k >= 0; k--) {
        tmp *= basis[gates[k]](x % 2, y % 2);
        x /= 2; y /= 2;
      }
      U(i, j) = tmp;
      U(i, j).reduce();
    }
  }
}

void Gate::tensor(Rmatrix & U, bool adj) const {
  if (!adj) this->tensor(U);
  else {
    int i, j, k, x, y;
    Elt tmp;

    for (i = 0; i < num_qubits; i++) {
      assert(gates[i] < basis_size + dim);
    }

    for (i = 0; i < dim; i++) {
      for (j = 0; j < dim; j++) {
        tmp = Elt(1, 0, 0, 0, 0);
        x = i; y = j;
        for (k = num_qubits-1; k >= 0; k--) {
          tmp *= basis[gates[k]](x % 2, y % 2);
          x /= 2; y /= 2;
        }
        U(j, i) = tmp.conj();
        U(j, i).reduce();
      }
    }
  }
}

void Gate::print() const {
  int i;
  char tmp;

  for (i = 0; i < num_qubits; i++) {
    tmp = gates[i];

    if (IS_C(tmp)) {
      cout << "C(" << (int)GET_TARGET(tmp) + 1 << ")";
    } else {
      assert(tmp <= basis_size);
      cout << " " << setw(3) << left << gate_names[tmp];
    }

    cout << "\n";
  }
}

/* Compute a unitary for the given gate */
void Gate::to_Rmatrix(Rmatrix & U, bool adj) const {
  int i;
  if (!config::tensors) {
    if (adj) {
      Gate G;
      this->adj(G);
      i = gate_hasher(G);
    } else {
      i = gate_hasher(*this);
    }
    if (i < (1 << 3*num_qubits) && gate_ht[i].rows() != 0) {
      U = gate_ht[i];
      return;
    }
  }

  for (i = 0; i < num_qubits; i++) {
    if (IS_C(gates[i])) {
      Gate A, B;
      Rmatrix V(dim, dim);
      char tmp = GET_TARGET(gates[i]);

      for (int j = 0; j < num_qubits; j++) {
        if (j == i) {
          A[j] = PROJ(i, 0);
          B[j] = PROJ(i, 1);
        } else if (j == tmp) {
          A[j] = I;
          B[j] = gates[j];
        } else {
          A[j] = gates[j];
          B[j] = gates[j];
        }
      }

      A.to_Rmatrix(V);
      B.to_Rmatrix(U);
      U += V;
      return;
    }
  }
  this->tensor(U, adj);
}

void Gate::to_Unitary(Unitary & U) const {
  Rmatrix tmp(dim, dim);
  this->to_Rmatrix(tmp);
  tmp.to_Unitary(U);
}

void Gate::permute(Gate & G, const char * perm) const {
  char gt;

  for(int i = 0; i < num_qubits; i++) {
    gt = gates[i];
    if (IS_C(gt)) {
      G[perm[i]] = C(perm[GET_TARGET(gt)]);
    } else {
      G[perm[i]] = gt;
    }
  }
}

void Gate::permute_adj(Gate & G, const char * perm) const {
  char gt;

  for (int i = 0; i < num_qubits; i++) {
    gt = gates[i];
    if (IS_C(gt)) {
      G[perm[i]] = C(perm[GET_TARGET(gt)]);
    } else {
      G[perm[i]] = adjoint[gt];
    }
  }
}

void Gate::permute(Gate & G, int i) const {
  const char * perm = from_lexi(i);
  this->permute(G, perm);
}

void Gate::permute_adj(Gate & G, int i) const {
  const char * perm = from_lexi(i);
  this->permute_adj(G, perm);
}

void Gate::output(ofstream & out) const {
  out.write(gates, num_qubits);
}
void Gate::input(ifstream & in) {
  if (gates == NULL) {
    gates = new char[num_qubits];
  }
  in.read(gates, num_qubits);
}

/* Determine if there is a nontrivial pair of gates that multiply to the identity */
bool nontrivial_id(const Gate & A, const Gate & B) {
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
    /*
    if (!IS_C(circ->G[i])) {
      if (G[i] == X && circ->G[i] == X) {
        ret = true;
        for (j = 0; j < num_qubits; j++) {
          ret = ret && !((G[j] == C(i)) xor (circ->G[j] == C(i)));
        }
      } else { 
        ret = ret || (G[i] == adjoint[circ->G[i]]);
      }
    }
    */
  }
  return ret || (bool)tmp;
}

const char subset[] = {0, 1, 6, 0, 0, 2, 3, 4, 5};
unsigned int gate_hasher(const Gate & R) {
  unsigned int ret = 0;
  int i;
  for (i = 0; i < num_qubits; i++) {
    if (IS_C(R[i])) {
      ret += 7 << (3*i);
    } else if (IS_PROJ(R[i])) {
      return 1 << (3*num_qubits);
    } else {
      ret += subset[R[i]] << (3*i);
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
    Gate G;
    Rmatrix R(dim, dim);
    Rmatrix * tmp = new Rmatrix[(1 << (3*num_qubits))];

    for (int j = 0; j < num_qubits; j++) {
      G[j] = I;
    }
    G.to_Rmatrix(R);
    tmp[gate_hasher(G)] = R;

    while(!((++G).eye())) {
      G.to_Rmatrix(R);
      tmp[gate_hasher(G)] = R;
    }

    gate_ht = tmp;
    config::tensors = false;
  }

  /*--------------------------- Initializing the identitites */
  init_identities();
}


void test_gate() {
  Gate A;
  A[0] = H;
  A[1] = I;
  A[2] = X;

  Gate B(A);
  assert(B == A);
  B[0] = I;
  assert(B != A);
  Gate * C = new Gate;
  (*C)[0] = I;
  (*C)[1] = I;
  (*C)[2] = I;
  assert(C->eye());
  A.adj(*C);
  const char * tmp = from_lexi(2);
  A.permute(B, tmp);

  Rmatrix R(dim, dim);
  Unitary U(dim, dim);
  B.to_Rmatrix(R);
  C->to_Unitary(U);

  delete C;
}

