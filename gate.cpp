#include "gate.h"

int num_qubits = 0;
int dim        = 0;
int num_swaps  = 0;
int num_weyl   = 0;

const string gate_names[] = {
  " I  ",
  " H  ",
  " X  ",  
  " Y  ",
  " Z  ",  
  " S  ",
  " S* ",
  " T  ",
  " T* ",
};

const char adjoint[] = {
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


Rmatrix * basis;
Rmatrix * swaps;
Unitary * weyl;
LaGenMatComplex subspace;

Unitary * maxU;

int fac(int n) {
  int ret = 1, i;
  for(i = n; i>1; i--) ret *= i;
  return ret;
}

int to_lexi(char * perm) {
  int i, tmp, ret = 0;
  for (i = 0 ; i < num_qubits-1; i++) {
    if (i == 0 || perm[i-1] <= perm[i]) {
      tmp = (perm[i] - i)*(num_qubits - 1 - i);
    } else {
      tmp = perm[i]*(num_qubits - 1 - i);
    }
    if (tmp < 0) tmp = 0;
    ret += tmp;
  }
  return ret;
}

char * from_lexi(int n) {
  char * ret = new char[num_qubits];
  int i, j, tmp, acc = n;
  for (i = 0; i < num_qubits; i++) {
    ret[i] = i;
  }
  for (i = 0; i < num_qubits-1; i++) {
    if (i == 0) {
      tmp = acc / (num_qubits - 1 - i);
    } else {
      tmp = acc / (num_qubits - 1 - i);
      if (tmp >= ret[i-1]) tmp += 1;
    }
    for (j = 0; j < num_qubits; j++) {
      if (ret[j] == tmp) {
        ret[j] = ret[i];
        ret[i] = tmp;
        break;
      }
    }
    acc = acc % (num_qubits - 1 - i);
  }
  return ret;
}


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
  for (i = 0; i < num_qubits; i++) {
    x = gates[i];
    if (IS_C(x)) {
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

/* --------------- Circuits */
Circuit::Circuit() { next = NULL; }
void delete_circuit(Circuit * circ) {
  Circuit * tmp;
  while (circ != NULL) {
    tmp = circ;
    circ = circ->next;
    delete tmp;
  }
}

void Circuit::print() const {
  const Circuit * tmp;
  char g;

  for (int i = 0; i < num_qubits; i++) {
    tmp = this;
    while (tmp) {
      g = tmp->G[i];
      if (IS_C(g)) {
        cout << "C(" << (int)GET_TARGET(g) + 1 << ")";
      } else {
        assert(g <= basis_size);
        cout << gate_names[g];
      }
      tmp = tmp->next;
    }
    cout << "\n";
  }
}

void Circuit::print(Circuit * snd) const {
  const Circuit * tmp;
  char g;

  for (int i = 0; i < num_qubits; i++) {
    tmp = this;
    while (tmp) {
      g = tmp->G[i];
      if (IS_C(g)) {
        cout << "C(" << (int)GET_TARGET(g) + 1 << ")";
      } else {
        assert(g <= basis_size);
        cout << gate_names[g];
      }
      tmp = tmp->next;
    }
    tmp = snd;
    while (tmp) {
      g = tmp->G[i];
      if (IS_C(g)) {
        cout << "C(" << (int)GET_TARGET(g) + 1 << ")";
      } else {
        assert(g <= basis_size);
        cout << gate_names[g];
      }
      tmp = tmp->next;
    }
    cout << "\n";
  }
}

Circuit * Circuit::adj(Circuit * last) const {
  Circuit * tmp = new Circuit;
  tmp->next = last;
  G.adj(tmp->G);

  if (next != NULL) return next->adj(tmp);
  else      return tmp;
}

Circuit * Circuit::permute(char * perm) const {
  Circuit * ret = new Circuit;
  G.permute(ret->G, perm);
  if (next == NULL) ret->next = NULL;
  else ret->next = next->permute(perm);
  return ret;
}

Circuit * Circuit::permute(int i) const {
  return this->permute(from_lexi(i));
}

void Circuit::to_Rmatrix(Rmatrix & U) const {
  if (next != NULL) {
	  Rmatrix V(dim, dim);
    G.to_Rmatrix(U);
    next->to_Rmatrix(V);
    U *= V;
  } else {
	  G.to_Rmatrix(U);
  }
}

void Circuit::to_Unitary(Unitary & U) const {
  Rmatrix tmp(dim, dim);
  this->to_Rmatrix(tmp);
  tmp.to_Unitary(U);
}

/*---------------------------------*/

void swap_qubits(char * perm, Rmatrix & swap) {
  int i, j, acc, index;
  char * col = new char[num_qubits];
  for (i = 0; i < dim; i++) {
    index = 0;
    acc = i;
    for (j = 0; j < num_qubits; j++) {
      col[j] = acc / pow(2, num_qubits - 1 - j);
      acc %= (int)pow(2, num_qubits - 1 - j);
    }

    for (j = 0; j < num_qubits; j++) {
      index += col[perm[j]] * (pow(2, num_qubits - 1 - j));
    }

    swap(index, i) = Elt(1, 0, 0, 0, 0);
  }
}

void init(int n) {
  num_qubits = n;
  num_swaps = fac(n);
  dim = pow(2, n);
  num_weyl = dim*dim;
  basis = new Rmatrix[basis_size + dim];
  swaps = new Rmatrix[num_swaps];
  weyl  = new Unitary[num_weyl];

  int i;

	for (i = 1; i < basis_size + dim; i++) {
	  basis[i] = zero(2, 2);
	}

  basis[I] = eye(2, 2);

	basis[H](0, 0) = Elt(0, 1, 0, -1, 1);
	basis[H](1, 0) = Elt(0, 1, 0, -1, 1);
	basis[H](0, 1) = Elt(0, 1, 0, -1, 1);
	basis[H](1, 1) = Elt(0, -1, 0, 1, 1);

	basis[X](0, 1) = Elt(1, 0, 0, 0, 0);
	basis[X](1, 0) = Elt(1, 0, 0, 0, 0);

	basis[Y](0, 1) = Elt(0, 0, -1, 0, 0);
	basis[Y](1, 0) = Elt(0, 0, 1, 0, 0);

	basis[Z](0, 0) = Elt(1, 0, 0, 0, 0);
	basis[Z](1, 1) = Elt(-1, 0, 0, 0, 0);

	basis[S](0, 0) = Elt(1, 0, 0, 0, 0);
	basis[S](1, 1) = Elt(0, 0, 1, 0, 0);

  basis[S].adj(basis[Sd]);

	basis[T](0, 0) = Elt(1, 0, 0, 0, 0);
	basis[T](1, 1) = Elt(0, 1, 0, 0, 0);

  basis[T].adj(basis[Td]);

  for (i = 0; i < dim; i++) {
    basis[PROJ((i / 2), (i % 2))](i%2, i%2) = Elt(1, 0, 0, 0, 0);
  }

  for (i = 0; i < num_swaps; i++) {
    swaps[i] = zero(dim, dim);
    swap_qubits(from_lexi(i), swaps[i]);
    /*
    t1 = from_lehmer(i);
    t2 = to_lehmer(t1);
    cout << t2 << " " << i << ": ";
    for (int j = 0; j < num_qubits; j++) {
      cout << (int)t1[j] << ",";
    }
    cout << "\n";
    */
  }

  Unitary x[dim];
  Unitary z[dim];
  int tmp, a, b, j;
  for (i = 0; i < dim; i++) {
    x[i] = Unitary::zeros(dim);
    z[i] = Unitary::zeros(dim);
    for (j = 0; j < dim; j++) {
      tmp = (2*PI*((i + j) % dim))/dim;
      (x[i])((i + j) % dim, j) = LaComplex(1, 0);
      (z[i])(j, j) = LaComplex(cos(tmp), sin(tmp));
    }
  }

  for (i = 0; i < num_weyl; i++) {
    a = i / dim;
    b = i % dim;
    weyl[i] = Unitary::zeros(dim);
    Blas_Mat_Mat_Mult(x[a], z[b], weyl[i], false, false, 1, 0);
  }

  /* Define the subspace */
  subspace = Unitary::rand(dim, SUBSPACE_SIZE);
  /*
  subspace = Unitary::zeros(dim, SUBSPACE_SIZE);
  subspace(0, 0) = LaComplex(1/sqrt(2), 0);
  subspace(7, 0) = LaComplex(1/sqrt(2), 0);
  subspace(1, 1) = LaComplex(1/2, 0);
  subspace(3, 1) = LaComplex(1/2, 0);
  subspace(4, 1) = LaComplex(-1/2, 0);
  subspace(6, 1) = LaComplex(-1/2, 0);
  */

  maxU = new Unitary(dim, dim);
  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      (*maxU)(i, j) = LaComplex(numeric_limits<double>::max(), numeric_limits<double>::max());
    }
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

double dist(const Rmatrix & M, const Rmatrix & N) {
  Unitary U(dim, dim);
  Unitary V(dim, dim);
  M.to_Unitary(U);
  N.to_Unitary(V);
  #ifndef PHASE
    return spec_norm(U - V);
  #else
    double acc = 0;
    Unitary A(dim, dim);
    Unitary B(dim, dim);

    for (int i = 0; i < num_weyl; i++) {
      Blas_Mat_Mat_Mult(weyl[i], U, A, false, true, 1, 0);
      Blas_Mat_Mat_Mult(U, A, A, false, false, 1, 0);
      Blas_Mat_Mat_Mult(weyl[i], V, B, false, true, 1, 0);
      Blas_Mat_Mat_Mult(V, B, B, false, false, 1, 0);
      acc += pow(spec_norm(A - B), 2);
    }

    return sqrt(acc);
  #endif
}
double dist(const Unitary & U, const Unitary & V) {
  #ifndef PHASE
    return spec_norm(U - V);
  #else
    double acc = 0;
    Unitary A(dim, dim);
    Unitary B(dim, dim);

    for (int i = 0; i < num_weyl; i++) {
      Blas_Mat_Mat_Mult(weyl[i], U, A, false, true, 1, 0);
      Blas_Mat_Mat_Mult(U, A, A, false, false, 1, 0);
      Blas_Mat_Mat_Mult(weyl[i], V, B, false, true, 1, 0);
      Blas_Mat_Mat_Mult(V, B, B, false, false, 1, 0);
      acc += pow(spec_norm(A - B), 2);
    }

    return sqrt(acc);
  #endif
}

/*
bool cmp_hash::operator()(const hash_t & a, const hash_t & b) {
  return a < b;
}

hash_t Hash_Unitary(const Unitary &U) {
  return (hash_t)(10000000*spec_norm(U - LaGenMatComplex::eye(dim)));
}

hash_t Hash_Rmatrix(const Rmatrix &U) {
  Unitary V = U.to_Unitary();
  return (hash_t)spec_norm(V - LaGenMatComplex::eye(dim));
}
*/

bool cmp_hash::operator()(const hash_t & a, const hash_t & b) {
  int i, j;
  double rea, reb, ima, imb;
  for (i = 0; i < SUBSPACE_SIZE; i++) {
    for (j = 0; j < SUBSPACE_SIZE; j++) {
      rea = real((LaComplex)a(i, j));
      ima = imag((LaComplex)a(i, j));
      reb = real((LaComplex)b(i, j));
      imb = imag((LaComplex)b(i, j));
      if (rea < reb || (rea == reb && ima < imb))
        return true;
      else if (rea > reb || ima > imb)
        return false;
    }
  }
        
  return false;
}

bool operator<(const hash_t & a, const hash_t & b) {
  struct cmp_hash x;
  return x(a, b);
}

bool operator==(const hash_t & a, const hash_t & b) {
  struct cmp_hash x;
  return !(x(a, b) || x(b, a));
}

hash_t Hash_Unitary(const Unitary & U) {
  int i, j;
  LaGenMatComplex tmp(dim, SUBSPACE_SIZE);
  Unitary V(SUBSPACE_SIZE, SUBSPACE_SIZE);

  Blas_Mat_Mat_Mult(U, subspace, tmp, false, false, 1, 0);
  Blas_Mat_Mat_Mult(subspace, tmp, V, true, false, 1, 0);
  for (i = 0; i < SUBSPACE_SIZE; i++) {
    for (j = 0; j < SUBSPACE_SIZE; j++) {
      V(i, j) = LaComplex(floor(PRECISION*real((LaComplex)V(i, j))), floor(PRECISION*imag((LaComplex)V(i, j))));
    }
  }

  return V;
}

hash_t Hash_Rmatrix(const Rmatrix & R) {
  int i, j;
  Unitary U = R.to_Unitary();
  LaGenMatComplex tmp(dim, SUBSPACE_SIZE);
  Unitary V(SUBSPACE_SIZE, SUBSPACE_SIZE);

  Blas_Mat_Mat_Mult(U, subspace, tmp, false, false, 1, 0);
  Blas_Mat_Mat_Mult(subspace, tmp, V, true, false, 1, 0);
  for (i = 0; i < SUBSPACE_SIZE; i++) {
    for (j = 0; j < SUBSPACE_SIZE; j++) {
      V(i, j) = LaComplex(floor(PRECISION*real((LaComplex)V(i, j))), floor(PRECISION*imag((LaComplex)V(i, j))));
    }
  }

  return V;
}

void permute(const Rmatrix & U, Rmatrix & V, int i) {
  V = U * swaps[i];
}

Canon canonicalize(const Rmatrix & U) {
  int i, j;
  hash_t d, min = *maxU;
  Rmatrix V(dim, dim), Vadj(dim, dim), best(dim, dim);
  Elt phase(0, 1, 0, 0, 0);

  Canon acc;

  for(i = 0; i < num_swaps; i++) {
    V = U;
    V *= swaps[i];
    V.adj(Vadj);

    for(j = 0; j < 8; j++) {
      if (j != 0) {
        V *= phase;
        Vadj *= phase;
      }
      d = Hash_Rmatrix(V);
      if (d < min) {
        min = d;
        best = V;
        acc.clear();
        acc.push_front({V, d, 0, i});
      } else if (!best.phase_eq(V) && d == min) {
        acc.push_front({V, d, 0, i});
      } 

      d = Hash_Rmatrix(Vadj);
      if (d < min) {
        min = d;
        best = Vadj;
        acc.clear();
        acc.push_front({Vadj, d, 1, i});
      } else if (!best.phase_eq(Vadj) && d == min) {
        acc.push_front({Vadj, d, 1, i});
      }
    }
  }

  return acc;
} 

void test() {
  init(1);

  Circuit * A = new Circuit;
  Circuit * B = new Circuit;
  Circuit * C = new Circuit;
  Circuit * D = new Circuit;
  Circuit * E = new Circuit;
  Circuit * F = new Circuit;
  Circuit * G = new Circuit;
  A->G[0] = H;
  B->G[0] = H;
  C->G[0] = H;
  D->G[0] = H;
  E->G[0] = H;
  F->G[0] = H;
  G->G[0] = H;

  A->next = B;
  B->next = C;
  C->next = D;
  D->next = E;
  E->next = F;
  F->next = G;
  G->next = NULL;

  Rmatrix U(dim, dim);
  Unitary V(dim, dim);
  A->to_Rmatrix(U);
  A->to_Unitary(V);
  U.print();
  cout << V << "\n";
}
