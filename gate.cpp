#include "gate.h"
#include <unordered_map>

unordered_map<Gate, Rmatrix, gate_hasher, gate_eq> gate_ht;
typedef pair<Gate, Rmatrix> gate_ht_elt;
bool use_ht = false;

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
    if (x == X) {
      flg = false;
      for(j = 0; j < num_qubits; j++) {
        flg = flg || (IS_C(gates[j]) && GET_TARGET(gates[j]) == i);
      }
      if (flg == false) return false;
    } else if (x == Y || x == Z) {
      return false;
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
      cout << gate_names[tmp];
    }

    cout << "\n";
  }
}

/* Compute a unitary for the given gate */
void Gate::to_Rmatrix(Rmatrix & U) const {
  int i;
  if (use_ht) {
    auto res = gate_ht.find(*this);
    if (res != gate_ht.end()) {
      U = res->second;
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
  this->tensor(U);
}

void Gate::to_Rmatrix(Rmatrix & U, bool adj) const {
  if (!adj) this->to_Rmatrix(U);
  else {
    int i;
    if (use_ht) {
      Gate G;
      this->adj(G);
      auto res = gate_ht.find(G);
      if (res != gate_ht.end()) {
        U = res->second;
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
}

void Gate::to_Unitary(Unitary & U) const {
  Rmatrix tmp(dim, dim);
  this->to_Rmatrix(tmp);
  tmp.to_Unitary(U);
}

void Gate::permute(Gate & G, char * perm) const {
  /* Do the permutation */
  for(int i = 0; i < num_qubits; i++) {
    G[i] = gates[perm[i]];
  }
  /* Fix the controls */
  for(int i = 0; i < num_qubits; i++) {
    if (IS_C(G[i])) {
      for (int j = 0; j < num_qubits; j++) {
        if (perm[j] == GET_TARGET(G[i])) {
          G[i] = C(j);
          break;
        }
      }
    }
  }
}

void Gate::permute_adj(Gate & G, char * perm) const {
  int i, j;

  for (i = 0; i < num_qubits; i++) {
    // Do permutation
    G[i] = gates[perm[i]];
    if (!IS_C(G[i])) {
      G[i] = adjoint[G[i]];
    }
  }
  // Fix controls
  for (i = 0; i < num_qubits; i++) {
    if (IS_C(G[i])) {
      for (j = 0; j < num_qubits; j++) {
        if (perm[j] == GET_TARGET(G[i])) {
          G[i] = C(j);
          break;
        }
      }
    }
  }
}

void Gate::permute(Gate & G, int i) const {
  char * perm = from_lexi(i);
  this->permute(G, perm);
  delete [] perm;
}

void Gate::permute_adj(Gate & G, int i) const {
  char * perm = from_lexi(i);
  this->permute_adj(G, perm);
  delete [] perm;
}

unsigned int gate_hasher::operator()(const Gate & R) const {
  unsigned int ret = 0;
  int i;
  for (i = 0; i < num_qubits; i++) {
    if (IS_C(R[i])) {
      ret += 9 * pow(10, i);
    } else {
      ret += R[i] * pow(10, i);
    }
  }
  return ret;
}

void init_ht() {
  Gate G;
  Rmatrix R(dim, dim);
  gate_ht.reserve(252);

  for (int j = 0; j < num_qubits; j++) {
    G[j] = I;
  }
  G.to_Rmatrix(R);
  gate_ht.insert(gate_ht_elt(G, R));

  while(!((++G).eye())) {
    G.to_Rmatrix(R);
    gate_ht.insert(gate_ht_elt(G, R));
  }

  use_ht = true;
}


void gate_test() {
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
  char * tmp = from_lexi(2);
  A.permute(B, tmp);
  delete [] tmp;

  Rmatrix R(dim, dim);
  Unitary U(dim, dim);
  B.to_Rmatrix(R);
  C->to_Unitary(U);

  delete C;
}

