#include "gate.h"

int num_qubits = 0;
int dim        = 0;
int reduced_dim= 0;
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

const int gatecost[2][9] = {{0,0,0,0,0,0,0,10,10}, 
                            {0,10,0,0,0,40,40,1000,1000}};
const int cnotcost[2] = {4,5};


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

/* --------------- Circuits */
Circuit::Circuit() { next = NULL; }
Circuit::Circuit(char g, Circuit * nxt) { G[0] = g; next = nxt; }
Circuit::Circuit(const Gate & Ga, Circuit * nxt) { G = Ga; next = nxt; }
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

Circuit * Circuit::reverse(Circuit * last) const {
  Circuit * ret = new Circuit;
  ret->next = last;
  ret->G = G;

  if (next != NULL) return next->adj(ret);
  else return ret;
}

Circuit * Circuit::adj(Circuit * last) const {
  Circuit * ret = new Circuit;
  ret->next = last;
  G.adj(ret->G);

  if (next != NULL) return next->adj(ret);
  else return ret;
}

Circuit * Circuit::permute(char * perm) const {
  Circuit * ret = new Circuit;
  G.permute(ret->G, perm);
  if (next == NULL) ret->next = NULL;
  else ret->next = next->permute(perm);
  return ret;
}

Circuit * Circuit::permute(int i) const {
  char * perm = from_lexi(i);
  Circuit * ret = this->permute(perm);
  delete [] perm;
  return ret;
}

Circuit * Circuit::permute_adj(char * perm, Circuit * last) const {
  Circuit * ret = new Circuit;
  ret->next = last;
  G.permute_adj(ret->G, perm);

  if (next != NULL) return next->permute_adj(perm, ret);
  else return ret;

  return ret;
}

Circuit * Circuit::permute_adj(int i, Circuit * last) const {
  char * perm = from_lexi(i);
  Circuit * ret = this->permute_adj(perm, last);
  delete [] perm;
  return ret;
}

Circuit * Circuit::append(Circuit * C) const {
  Circuit * ret = new Circuit;
  ret->G = G;
  if (next == NULL) {
    ret->next = C;
  } else {
    ret->next = next->append(C);
  }
  return ret;
}

const Gate & Circuit::last() const {
  const Circuit * tmp = this;
  while (tmp->next != NULL) tmp = tmp->next;
  return tmp->G;
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

void Circuit::to_Rmatrix(Rmatrix & U, bool adj) const {
  if (!adj) this->to_Rmatrix(U);
  else {
    if (next != NULL) {
      Rmatrix V(dim, dim);
      G.to_Rmatrix(V);
      next->to_Rmatrix(U, adj);
      U.left_multiply(V);
    } else {
      G.to_Rmatrix(U);
    }
  }
}

void Circuit::to_Unitary(Unitary & U) const {
  Rmatrix tmp(dim, dim);
  this->to_Rmatrix(tmp);
  tmp.to_Unitary(U);
}

int Circuit::cost(Arch a) const {
  const Circuit * tmp = this;
  int ret = 0, i;

  switch (a) {
    case STEANE:
      break;
    case SURFACE:
      while (this != NULL) {
        for (i = 0; i < num_qubits; i++) {
          if (IS_C(tmp->G[i])) {
            ret += cnotcost[a];
          } else {
            ret += gatecost[a][tmp->G[i]];
          }
        }
        tmp = tmp->next;
      }
      break;
    default:
      break;
  }
  return ret;
}

/*---------------------------------*/

int swap_i(char * perm, int i) {
  char * x = new char[num_qubits], * y = new char[num_qubits];
  int j, ret = i;

  /* Determine the basis element corresponding to i */
  for (j = 0; j < num_qubits; j++) {
    x[j] = ret / (int)pow(2, num_qubits - 1 - j);
    ret %= (int)pow(2, num_qubits - 1 - j);
  }

  /* Determine the permuted basis element */
  for (j = 0; j < num_qubits; j++) {
    y[perm[j]] = x[j];
  }

  ret = 0;
  /* Compute the integer representing the permuted basis element */
  for (j = num_qubits - 1; j >= 0; j--) {
    ret += y[j] * (pow(2, num_qubits - 1 - j));
  }

  delete [] x;
  delete [] y;

  return ret;
}

void swap_qubits(char * perm, Rmatrix & swap) {
  int i, j, ip, jp;
  for (i = 0; i < dim; i++) {
    ip = swap_i(perm, i);
    swap(ip, i) = Elt(1, 0, 0, 0, 0);
  }
}

void init(int n, int m) {
  num_qubits = n;
  num_swaps = fac(n);
  dim = pow(2, n);
  reduced_dim = pow(2, m);
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

  char * swap_tmp;
  for (i = 0; i < num_swaps; i++) {
    swaps[i] = zero(dim, dim);
    swap_tmp = from_lexi(i);
    swap_qubits(swap_tmp, swaps[i]);
    delete [] swap_tmp;
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
  int i, j, m = R.rows() - 1, n = R.cols() - 1;
  Unitary U = R.to_Unitary();
  LaGenMatComplex tmp(R.rows(), SUBSPACE_SIZE);
  Unitary V(SUBSPACE_SIZE, SUBSPACE_SIZE);

  Blas_Mat_Mat_Mult(U, subspace(LaIndex(0, n), LaIndex(0, SUBSPACE_SIZE-1)), tmp, false, false, 1, 0);
  Blas_Mat_Mat_Mult(subspace(LaIndex(0, m), LaIndex(0, SUBSPACE_SIZE-1)), tmp, V, true, false, 1, 0);
  for (i = 0; i < SUBSPACE_SIZE; i++) {
    for (j = 0; j < SUBSPACE_SIZE; j++) {
      V(i, j) = LaComplex(PRECISION*real((LaComplex)V(i, j)), PRECISION*imag((LaComplex)V(i, j)));
    }
  }

  return V;
} 

void permute(const Rmatrix & U, Rmatrix & V, int i) {
  Rmatrix tmp(dim, dim);
  swaps[i].adj(tmp);
  V = tmp * U * swaps[i];
}

/* Returns the canonical form(s), ie. the lowest hashing unitary for
   each permutation, inversion, and phase factor */
Canon canonicalize(const Rmatrix & U, bool sym) {
  int i, j;
  hash_t d, min = *maxU;
  Rmatrix V(dim, dim), Vadj(dim, dim), best(dim, dim);
  Elt phase(0, 1, 0, 0, 0);

  Canon acc;

  if (sym) {
    for (i = 0; i < num_swaps; i++) {
      permute(U, V, i);
      V.adj(Vadj);

      for (j = 0; j < 8; j++) {
        if (j != 0) {
          V *= phase;
          Vadj *= phase;
        }

        d = Hash_Rmatrix(V);
        if (d < min) {
          min = d;
          best = V;
          acc.clear();
          acc.push_front({V, d, false, i});
        } else if (!best.phase_eq(V) && d == min) {
          acc.push_front({V, d, false, i});
        } 

        d = Hash_Rmatrix(Vadj);
        if (d < min) {
          min = d;
          best = Vadj;
          acc.clear();
          acc.push_front({Vadj, d, true, i});
        } else if (!best.phase_eq(Vadj) && d == min) {
          acc.push_front({Vadj, d, true, i});
        }
      }
    }
  } else {
#ifdef PHASE
    V = U;
    for (j = 0; j < 8; j++) {
      if (j != 0) {
        V *= phase;
      }

      d = Hash_Rmatrix(V);
      if (d < min) {
        min = d;
        best = V;
        acc.clear();
        acc.push_front({V, d, false, i});
      } else if (!best.phase_eq(V) && d == min) {
        acc.push_front({V, d, false, i});
      }
    }
#else 
    acc.push_front({U, Hash_Rmatrix(U), false, 0});
#endif
  }

  return acc;
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

void circuit_test() {
  Circuit * A = new Circuit;
  Circuit * B = new Circuit;
  B->next = A;

  A->G[0] = H;
  A->G[1] = C(2);
  A->G[2] = T;
  B->G[0] = X;
  B->G[1] = Y;
  B->G[2] = Z;

  Circuit * C = B->reverse(NULL);
  Circuit * D = B->adj(NULL);
  Circuit * E = B->permute(2);
  Circuit * F = C->append(D);

  
  Rmatrix a(dim, dim), b(dim, dim), c(dim, dim);
  B->to_Rmatrix(a);
  D->to_Rmatrix(b);
  a.adj(c);
  assert(c == b);
  E->to_Rmatrix(b);
  permute(a, c, 2);
  assert(c == b);

  Canon canon = canonicalize(a, true);

  delete A;
  delete B;
  delete_circuit(C);
  delete_circuit(E);
  delete_circuit(F);
}
