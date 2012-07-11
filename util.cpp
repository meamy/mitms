#include "util.h"
#include <pthread.h>

pthread_mutex_t blas_lock;

const Unitary * weyl;
const subs_t * subspace;
const subs_t * subspace_adj;
const hash_t * maxU;

/*------------------- Utilities */
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

/*------------------- Distance calculations & norms */
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
  if (config::mod_phase) {
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
  } else {
    return spec_norm(U - V);
  }
}
double dist(const Unitary & U, const Unitary & V) {
  if (config::mod_phase) {
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
  } else {
    return spec_norm(U - V);
  }
}

/*--------------------- Keys & hashing utilities */
bool cmp_hash::operator()(const hash_t & a, const hash_t & b) const {
  int i, j;
  int rsgn, isgn;
  double rea, reb, ima, imb;
  for (i = 0; i < config::key_dimension; i++) {
    for (j = 0; j < config::key_dimension; j++) {
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

bool eq_hash::operator()(const hash_t & a, const hash_t & b) const {
  int i, j;
  double rea, reb, ima, imb;
  for (i = 0; i < config::key_dimension; i++) {
    for (j = 0; j < config::key_dimension; j++) {
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

bool operator<(const hash_t & a, const hash_t & b) {
  struct cmp_hash x;
  return x(a, b);
}

bool operator==(const hash_t & a, const hash_t & b) {
  struct eq_hash x;
  return x(a, b);
}

// Hash function for unordered maps 
unsigned int hasher::operator()(const hash_t & a) const {
  int sz = config::key_dimension*config::key_dimension;
  int x, y, z;
  int tmp, i, j;
  unsigned int ret = 0;

  for (i = 0; i < sz; i++) {
    x = i / config::key_dimension;
    y = i % config::key_dimension;
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
  LaGenMatComplex tmp(dim, config::key_dimension);
  hash_t V(config::key_dimension, config::key_dimension);

  pthread_mutex_lock(&blas_lock);
  Blas_Mat_Mat_Mult(U, *subspace, tmp, false, false, 1, 0);
  Blas_Mat_Mat_Mult(*subspace, tmp, V, true, false, 1, 0);
  pthread_mutex_unlock(&blas_lock);

  if (config::approximate) {
    for (i = 0; i < config::key_dimension; i++) {
      for (j = 0; j < config::key_dimension; j++) {
        V(i, j) = LaComplex(
            floor(config::precision*real((LaComplex)V(i, j))), 
            floor(config::precision*imag((LaComplex)V(i, j))));
      }
    }
  }

  return V;
}

hash_t Hash_Rmatrix(const Rmatrix & R) {
  int i, j, m = R.rows() - 1, n = R.cols() - 1;
  const Unitary U = config::mod_phase ? R.to_Unitary_canon() : R.to_Unitary();
  LaGenMatComplex tmp(R.rows(), config::key_dimension);
  hash_t V(config::key_dimension, config::key_dimension);

  pthread_mutex_lock(&blas_lock);
  Blas_Mat_Mat_Mult(U, (*subspace)(LaIndex(0, n), LaIndex(0, config::key_dimension-1)), tmp, false, false, 1, 0);
  Blas_Mat_Mat_Mult((*subspace)(LaIndex(0, m), LaIndex(0, config::key_dimension-1)), tmp, V, true, false, 1, 0);
  pthread_mutex_unlock(&blas_lock);

  if (config::approximate) {
    for (i = 0; i < config::key_dimension; i++) {
      for (j = 0; j < config::key_dimension; j++) {
        V(i, j) = LaComplex(
            floor(config::precision*real((LaComplex)V(i, j))), 
            floor(config::precision*imag((LaComplex)V(i, j))));
      }
    }
  }

  return V;
}

/* Returns the canonical form(s), ie. the lowest hashing unitary for
   each permutation, inversion, and phase factor */
Canon * canonicalize(const Rmatrix & U, bool phse, bool perms, bool invs) {
  int i, j;
  int ph = 1;
  int pe = perms ? num_perms : 1;
  hash_t d, min = *maxU;
  Rmatrix V(dim, dim), Vadj(dim, dim), best(dim, dim);
  Elt phase(0, 1, 0, 0, 0);

  Canon * acc = new Canon;
  struct triple ins;

  for (i = 0; i < pe; i++) {
    U.permute(V, i);
    V.adj(Vadj);

    for (j = 0; j < ph; j++) {
      if (j != 0) {
        V *= phase;
        Vadj *= phase;
      }

      d = Hash_Rmatrix(V);
      if (d < min) {
        min = d;
        best = V;
        acc->clear();
        acc->push_front(triple(V, d, false, i));
      } else if (!best.phase_eq(V) && d == min) {
        acc->push_front(triple(V, d, false, i));
      } 

      if (invs) {
        d = Hash_Rmatrix(Vadj);
        if (d < min) {
          min = d;
          best = Vadj;
          acc->clear();
          acc->push_front(triple(Vadj, d, true, i));
        } else if (!best.phase_eq(Vadj) && d == min) {
          acc->push_front(triple(Vadj, d, true, i));
        }
      }
    }
  }

  return acc;
}

void output_key(ofstream & out, const hash_t & key) {
  int i, j;
  double * c = new double[2];
  for (i = 0; i < config::key_dimension; i++) {
    for (j = 0; j < config::key_dimension; j++) {
      c[0] = real((LaComplex)key(i, j));
      c[1] = imag((LaComplex)key(i, j));
      out.write((char *)c, 2*sizeof(double));
    }
  }
  delete [] c;
}
          
void input_key (ifstream & in, hash_t & key) {
  int i, j;
  double * c = new double[2];
  if ((key.rows() | key.cols()) != config::key_dimension) {
    key = hash_t(config::key_dimension, config::key_dimension);
  }
  for (i = 0; i < config::key_dimension; i++) {
    for (j = 0; j < config::key_dimension; j++) {
      in.read((char *)c, 2*sizeof(double));
      key(i, j) = LaComplex(c[0], c[1]);
    }
  }
  delete [] c;
}

/*---------------------------------*/

void init_util() {
  int i, j, a, b;
  int aa = numeric_limits<int>::max();
  double tmp;

  pthread_mutex_init(&blas_lock, NULL);

  /*----------------------- Initializing generalized Paulis */
  Unitary * weyl_tmp = new Unitary[num_weyl];

  Unitary x[dim];
  Unitary z[dim];
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
    weyl_tmp[i] = Unitary::zeros(dim);
    Blas_Mat_Mat_Mult(x[a], z[b], weyl_tmp[i], false, false, 1, 0);
  }
  weyl = weyl_tmp;

  /*----------------------- Initializing key subspace */
  subs_t * subspace_tmp = new subs_t;
  *subspace_tmp = subs_t::rand(dim, config::key_dimension);
  subspace = subspace_tmp;

  Unitary * maxU_tmp = new hash_t(dim, dim);
  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      (*maxU_tmp)(i, j) = LaComplex(numeric_limits<double>::max(), numeric_limits<double>::max());
    }
  }
  maxU = maxU_tmp;

  srand(time(NULL));
}

