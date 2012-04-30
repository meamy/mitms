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
Rmatrix * key_permutations;
Unitary * weyl;
subs_t subspace;
subs_t subspace_adj;

hash_t * maxU;

int sgn(double a, double b) {
  return (b < a) - (a < b);
}

int max (int a, int b) {
  if (a > b) return a;
  else return b;
}

bool double_eq(double a, double b) {
  double epsilon = numeric_limits<double>::epsilon();
  return a == b;
//  return abs(a - b) < 0.00000001;
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
      double s = spec_norm(A - B);
      acc += s*s;
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
      double s = spec_norm(A - B);
      acc += s*s;
    }

    return sqrt(acc);
  #endif
}

#if SUBSPACE_ABS
bool cmp_hash::operator()(const hash_t & a, const hash_t & b) const {
  return a < b;
}

bool eq_hash::operator()(const hash_t & a, const hash_t & b) const {
  return a == b;
}

unsigned int hasher::operator()(const hash_t & a) const {
  int sz = SUBSPACE_SIZE*SUBSPACE_SIZE;
  int x, y, z;
  int tmp, i, j;
  unsigned int ret = 0;

  for (i = 0; i < sz; i++) {
    x = i / SUBSPACE_SIZE;
    y = i % SUBSPACE_SIZE;
    z = (int)ceil((float)(9-i) / (float)sz);
    for (j = 0; j < z; j++) {
      const complex<double> w = a(x, y).to_complex();
      tmp = (int)((real(w)+imag(w))*pow(10.0, j+1)) % 10;
      ret += tmp*pow(10.0, 9-i*j);
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
 // hash_t U(SUBSPACE_SIZE, SUBSPACE_SIZE);
  return subspace_adj * R * subspace;

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
 // return U;
}

#else
bool cmp_hash::operator()(const hash_t & a, const hash_t & b) const {
  int i, j;
  int rsgn, isgn;
  double rea, reb, ima, imb;
  for (i = 0; i < SUBSPACE_SIZE; i++) {
    for (j = 0; j < SUBSPACE_SIZE; j++) {
      rea = real((LaComplex)a(i, j));
      reb = real((LaComplex)b(i, j));
      rsgn = sgn(rea, reb);
      if (rsgn == 0) {
        ima = imag((LaComplex)a(i, j));
        imb = imag((LaComplex)b(i, j));
        isgn = sgn(ima, imb);
        if (isgn != 0) {
          return (isgn < 0);
        }
      } else {
        return (rsgn < 0);
      }
      /*
      if (rea < reb || (rea == reb && ima < imb))
        return true;
      else if (rea > reb || ima > imb)
        return false;
        */
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
  int sz = SUBSPACE_SIZE*SUBSPACE_SIZE;
  int x, y, z;
  int tmp, i, j;
  unsigned int ret = 0;

  for (i = 0; i < sz; i++) {
    x = i / SUBSPACE_SIZE;
    y = i % SUBSPACE_SIZE;
    z = (int)ceil((float)(9-i) / (float)sz);
    for (j = 0; j < z; j++) {
      tmp = (int)((real((LaComplex)(a(x, y))) + imag((LaComplex)(a(x, y)))) * pow(10.0, j+1)) % 10;
      ret += tmp*pow(10.0, 9-i*j);
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
  /*
  for (i = 0; i < SUBSPACE_SIZE; i++) {
    for (j = 0; j < SUBSPACE_SIZE; j++) {
      V(i, j) = LaComplex(PRECISION*real((LaComplex)V(i, j)), PRECISION*imag((LaComplex)V(i, j)));
    }
  }
*/

  return V;
}
#endif

void permute_hash(const Rmatrix & U, Rmatrix & V, int i) {
  V = key_permutations[i + num_swaps] * U * key_permutations[i];
}

void permute_adj_hash(const Rmatrix & U, Rmatrix & V, int i) {
  Rmatrix tmp(dim, dim);
  tmp = key_permutations[i + num_swaps] * U * key_permutations[i];
  tmp.adj(V);
}

bool equiv(const Rmatrix & M, const Rmatrix & N) {
  int i;
  Rmatrix tmp(dim, dim), t(dim, dim);
  for (i = 0; i < 2*num_swaps; i++) {
    if (i % 2 == 0) {
      M.permute(t, i / 2);
    } else {
      M.permute_adj(t, i / 2);
    }
    if (t.phase_eq(N)) return true;
  }
  return false;
}

/* Returns the canonical form(s), ie. the lowest hashing unitary for
   each permutation, inversion, and phase factor */
#if SUBSPACE_ABS
Canon canonicalize(const hash_t & U, bool sym) {
  int i, j;
  hash_t d, min = *maxU, tmp;
  Rmatrix V, Vadj;
  Elt phase(0, 1, 0, 0, 0);

  Canon acc;

  if (sym) {
    for (i = 0; i < num_swaps; i++) {
      V = key_permutations[i + num_swaps] * U * key_permutations[i];
      V.adj(Vadj);

#if PHASE
      for (j = 0; j < 8; j++) {
        if (j != 0) {
          V *= phase;
          Vadj *= phase;
        }
#endif

        if (V < min) {
          min = V;
          acc.clear();
          acc.push_front({V, V, false, i});
        } else if (!min.phase_eq(V)) {
          acc.push_front({V, V, false, i});
        } 

        if (Vadj < min) {
          min = Vadj;
          acc.clear();
          acc.push_front({Vadj, Vadj, true, i});
        } else if (!min.phase_eq(Vadj)) {
          acc.push_front({Vadj, Vadj, true, i});
        }
#if PHASE
      }
#endif
    }
  } else {
    acc.push_front({U, U, false, 0});
  }

  return acc;
}
#else
Canon canonicalize(const Rmatrix & U, bool sym) {
  int i, j;
  hash_t d, min = *maxU, tmp;
  Rmatrix V(dim, dim), Vadj(dim, dim), best(dim, dim);
  Elt phase(0, 1, 0, 0, 0);

  Canon acc;
  struct triple ins;

  if (sym) {
    for (i = 0; i < num_swaps; i++) {
      U.permute(V, i);
      V.adj(Vadj);

#if PHASE
      for (j = 0; j < 8; j++) {
        if (j != 0) {
          V *= phase;
          Vadj *= phase;
        }
#endif

        d = Hash_Rmatrix(V);
        if (d < min) {
          min = d;
          best = V;
          acc.clear();
          acc.push_front(triple(V, d, false, i));
        } else if (!best.phase_eq(V) && d == min) {
          acc.push_front(triple(V, d, false, i));
        } 

        d = Hash_Rmatrix(Vadj);
        if (d < min) {
          min = d;
          best = Vadj;
          acc.clear();
          acc.push_front(triple(Vadj, d, true, i));
        } else if (!best.phase_eq(Vadj) && d == min) {
          acc.push_front(triple(Vadj, d, true, i));
        }
#if PHASE
      }
#endif
    }
  } else {
    acc.push_front(triple(U, Hash_Rmatrix(U), false, 0));
  }

  return acc;
}
#endif

/*---------------------------------*/

void init(int n, int m) {
  num_qubits = n;
  num_swaps = fac(n);
  dim = 1 << n;
  reduced_dim = 1 << m;
  num_weyl = dim*dim;

  int i;
  init_permutations(n);

  /*-------------------------- Initializing matrices */
  basis = new Rmatrix[basis_size + dim];
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

  /*----------------------- Initializing generalized Paulis */
  weyl  = new Unitary[num_weyl];

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

  /*----------------------- Initializing key subspace */
  subspace = subs_t::rand(dim, SUBSPACE_SIZE);
#if SUBSPACE_ABS
  subspace.adj(subspace_adj);
#endif

  maxU = new hash_t(dim, dim);
  int aa = numeric_limits<int>::max();
  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
#if SUBSPACE_ABS
      (*maxU)(i, j) = Elt(aa, aa, aa, aa, aa);
#else
      (*maxU)(i, j) = LaComplex(numeric_limits<double>::max(), numeric_limits<double>::max());
#endif
    }
  }

  /*----------------------- Initializing permutation matrices */
#if SUBSPACE_ABS
  key_permutations = new Rmatrix[2*num_swaps];
  for (i = 0; i < 2*num_swaps; i++) {
    key_permutations[i] = subspace_adj * col_permutation(dim, dim, i) * subspace;
  }
#endif

  srand(time(NULL));
}

