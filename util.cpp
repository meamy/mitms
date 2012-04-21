#include "util.h"

int num_qubits = 0;
int dim        = 0;
int reduced_dim= 0;
int num_swaps  = 0;
int num_weyl   = 0;

const string gate_names[basis_size] = {
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

const char adjoint[basis_size] = {
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

const int gate_cost[2][basis_size] = {{0,0,0,0,0,0,0,10,10}, 
                            {0,10,0,0,0,40,40,1000,1000}};
const int cnot_cost[2] = {4,5};

Rmatrix * basis;
Rmatrix * swaps;
Unitary * weyl;
subs_t subspace;
subs_t subspace_adj;

hash_t * maxU;

int max (int a, int b) {
  if (a > b) return a;
  else return b;
}

int fac(int n) {
  int ret = 1, i;
  for(i = n; i>1; i--) ret *= i;
  return ret;
}

bool double_eq(double a, double b) {
  double epsilon = numeric_limits<double>::epsilon();
  return a == b;
//  return abs(a - b) < 0.00000001;
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

#if SUBSPACE_ABS
bool cmp_hash::operator()(const hash_t & a, const hash_t & b) const {
  int i, j;
  for (i = 0; i < SUBSPACE_SIZE; i++) {
    for (j = 0; j < SUBSPACE_SIZE; j++) {
      if (a(i,j) < b(i,j)) return true;
      if (!double_eq(a(i, j), b(i, j))) return false;
    }
  }
  return false;
}

bool operator<(const hash_t & a, const hash_t & b) {
  struct cmp_hash x;
  return x(a, b);
}

bool eq_hash::operator()(const hash_t & a, const hash_t & b) const {
  int i, j;
  for (i = 0; i < SUBSPACE_SIZE; i++) {
    for (j = 0; j < SUBSPACE_SIZE; j++) {
      if (!double_eq(a(i, j), b(i, j))) return false;
    }
  }
  return true;
}

bool operator==(const hash_t & a, const hash_t & b) {
  struct eq_hash x;
  return x(a, b);
}

unsigned int hasher::operator()(const hash_t & a) const {
  int sz = pow(SUBSPACE_SIZE, 2);
  int x, y, z;
  int tmp, i, j;
  unsigned int ret = 0;

  for (i = 0; i < sz; i++) {
    x = i / SUBSPACE_SIZE;
    y = i % SUBSPACE_SIZE;
    z = (int)ceil((float)(9-i) / (float)sz);
    for (j = 0; j < z; j++) {
      tmp = (int)((a(x, y)) * pow(10, j+1)) % 10;
      ret += tmp*pow(10, 9-i*j);
    }
  }
  return ret;
}

hash_t Hash_Unitary(const Unitary & U) {
  hash_t W(dim, dim);
  /*
  int i, j;
  LaGenMatDouble W(dim, dim);
  for (i = 0; i < SUBSPACE_SIZE; i++) {
    for (j = 0; j < SUBSPACE_SIZE; j++) {
      W(i, j) = abs((LaComplex)(V(i, j)));
    }
  }
  LaGenMatDouble tmp(dim, SUBSPACE_SIZE);
  LaGenMatDouble V(SUBSPACE_SIZE, SUBSPACE_SIZE);
  hash_t W(SUBSPACE_SIZE, SUBSPACE_SIZE);

  Blas_Mat_Mat_Mult(U, subspace, tmp, false, false, 1, 0);
  Blas_Mat_Mat_Mult(subspace, tmp, V, true, false, 1, 0);
*/
  return W;
}

hash_t Hash_Rmatrix(const Rmatrix & R) {
  int i, j, m = R.rows() - 1, n = R.cols() - 1;
  hash_t U(SUBSPACE_SIZE, SUBSPACE_SIZE);

  (subspace_adj * R * subspace).to_Unitary_abs(U);

/*
  subs_t U(R.rows(), R.cols());
  R.to_Unitary_abs(U);
  subs_t tmp(R.rows(), SUBSPACE_SIZE);
  subs_t V(SUBSPACE_SIZE, SUBSPACE_SIZE);
  hash_t W(SUBSPACE_SIZE, SUBSPACE_SIZE);

  Blas_Mat_Mat_Mult(U, subspace(LaIndex(0, n), LaIndex(0, SUBSPACE_SIZE-1)), tmp, false, false, 1, 0);
  Blas_Mat_Mat_Mult(subspace(LaIndex(0, m), LaIndex(0, SUBSPACE_SIZE-1)), tmp, V, true, false, 1, 0);
  
  for (i = 0; i < SUBSPACE_SIZE; i++) {
    for (j = 0; j < SUBSPACE_SIZE; j++) {
      W(i, j) = abs((LaComplex)(V(i, j)));
    }
  }
*/
  return U;
}

#else
bool cmp_hash::operator()(const hash_t & a, const hash_t & b) const {
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

bool eq_hash::operator()(const hash_t & a, const hash_t & b) const {
  int i, j;
  double rea, reb, ima, imb;
  for (i = 0; i < SUBSPACE_SIZE; i++) {
    for (j = 0; j < SUBSPACE_SIZE; j++) {
      rea = real((LaComplex)a(i, j));
      ima = imag((LaComplex)a(i, j));
      reb = real((LaComplex)b(i, j));
      imb = imag((LaComplex)b(i, j));
      if (rea != reb || ima != imb)
        return false;
    }
  }
  return true;
}

bool operator==(const hash_t & a, const hash_t & b) {
  struct eq_hash x;
  return x(a, b);
}

unsigned int hasher::operator()(const hash_t & a) const {
  int sz = pow(SUBSPACE_SIZE, 2);
  int x, y, z;
  int tmp, i, j;
  unsigned int ret = 0;

  for (i = 0; i < sz; i++) {
    x = i / SUBSPACE_SIZE;
    y = i % SUBSPACE_SIZE;
    z = (int)ceil((float)(9-i) / (float)sz);
    for (j = 0; j < z; j++) {
      tmp = (int)((real((LaComplex)(a(x, y))) + imag((LaComplex)(a(x, y)))) * pow(10, j+1)) % 10;
      ret += tmp*pow(10, 9-i*j);
    }
  }
  return ret;
}

hash_t Hash_Unitary(const Unitary & U) {
  int i, j;
  LaGenMatComplex tmp(dim, SUBSPACE_SIZE);
  hash_t V(SUBSPACE_SIZE, SUBSPACE_SIZE);

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
  Unitary U = R.to_Unitary_canon();
  LaGenMatComplex tmp(R.rows(), SUBSPACE_SIZE);
  hash_t V(SUBSPACE_SIZE, SUBSPACE_SIZE);

  Blas_Mat_Mat_Mult(U, subspace(LaIndex(0, n), LaIndex(0, SUBSPACE_SIZE-1)), tmp, false, false, 1, 0);
  Blas_Mat_Mat_Mult(subspace(LaIndex(0, m), LaIndex(0, SUBSPACE_SIZE-1)), tmp, V, true, false, 1, 0);
  for (i = 0; i < SUBSPACE_SIZE; i++) {
    for (j = 0; j < SUBSPACE_SIZE; j++) {
      V(i, j) = LaComplex(PRECISION*real((LaComplex)V(i, j)), PRECISION*imag((LaComplex)V(i, j)));
    }
  }

 // cout << hasher(V) << "\n";

  return V;
}
#endif

void permute(const Rmatrix & U, Rmatrix & V, int i) {
  Rmatrix tmp(dim, dim);
  V = swaps[i + num_swaps] * U * swaps[i];
}

void permute_inv(const Rmatrix & U, Rmatrix & V, int i) {
  Rmatrix tmp(dim, dim);
  V = swaps[i] * U * swaps[i + num_swaps];
}

bool equiv(const Rmatrix & M, const Rmatrix & N) {
  int i;
  Rmatrix tmp(dim, dim), t(dim, dim);
  for (i = 0; i < 2*num_swaps; i++) {
    if (i % 2 == 0) {
      permute(M, t, i / 2);
    } else {
      M.adj(tmp);
      permute(tmp, t, i / 2);
    }
    if (t.phase_eq(N)) return true;
  }
  return false;
}

/* Returns the canonical form(s), ie. the lowest hashing unitary for
   each permutation, inversion, and phase factor */
Canon canonicalize(const Rmatrix & U, bool sym) {
  int i, j;
  hash_t d, min = *maxU, tmp;
  Rmatrix V(dim, dim), Vadj(dim, dim), best(dim, dim);
  Elt phase(0, 1, 0, 0, 0);

  Canon acc;

  if (sym) {
    for (i = 0; i < num_swaps; i++) {
      permute(U, V, i);
      V.adj(Vadj);

#if PHASE
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
#else
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
#endif
    }
  } else {
  /*  
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
   */ 
    acc.push_front({U, Hash_Rmatrix(U), false, 0});
  }

  return acc;
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
  swaps = new Rmatrix[2*num_swaps];
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
    swaps[i].adj(swaps[num_swaps + i]);
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
  subspace = subs_t::rand(dim, SUBSPACE_SIZE);
#if SUBSPACE_ABS
  subspace.adj(subspace_adj);
#endif

  maxU = new hash_t(dim, dim);
  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
#if SUBSPACE_ABS
      (*maxU)(i, j) = numeric_limits<double>::max();
#else
      (*maxU)(i, j) = LaComplex(numeric_limits<double>::max(), numeric_limits<double>::max());
#endif
    }
  }
  srand(time(NULL));
}

