#ifndef UTIL
#define UTIL

#define LA_COMPLEX_SUPPORT

#include "matrix.h"
#include <list>
#include <iostream>
#include <string>
#include <assert.h>
#include <map>
#include <limits>
#include <unordered_map>

#include <gmc.h>
#include <blas3pp.h>
#include <blas2pp.h>
#include <blaspp.h>
#include <laslv.h>
#include <lavc.h>

#define PI 3.14159

/* Number of basis gates */
#define basis_size 9

/* 1 qubit basis gates */
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
#define C(x)          (x | 0x40)
#define IS_C(x)       (x & 0x40)
#define GET_TARGET(x) (x & 0x3F)

/* Configs */
#define SUBSPACE_SIZE 1   // Size of subspace we store
#define SUBSPACE_ABS  false// take the absolute value of the subspace matrix
#define PRECISION 1       // Precision
#define PHASE false       // Whether we mod out phase
#define MAX_SEQ 50        // Max circuit length
#define CLIFF 50          // Max depth we try to find clifford circuits
#define SYMMS true        // Whether we mod out symmetries
#define CHECK_EQUIV false  // Whether we check to make sure two circuits are equiv
#define ORDERED true      // Whether we should use an ordered map
#define TENSORS true      // Whether to store gates as tensor products of gates
#define TDEPTH  false      // Whether we want to search by T-depth


using namespace std;

enum Arch { STEANE = 0, SURFACE = 1 };

extern int num_qubits;
extern int num_swaps;
extern int dim;
extern int reduced_dim;
extern const string gate_names[basis_size];
extern const char   adjoint[basis_size];
extern const int    gate_cost[][basis_size];
extern const int    cnot_cost[];

extern Rmatrix * basis;

typedef LaGenMatComplex Unitary;
typedef list< struct triple > Canon;
#if SUBSPACE_ABS
  typedef Rmatrix hash_t;
  typedef Rmatrix subs_t;
#else
  typedef LaGenMatComplex hash_t;
  typedef LaGenMatComplex subs_t;
#endif

struct triple {
  Rmatrix mat;
  hash_t key;
  bool adjoint;
  int permutation;
};

int max (int a, int b);
char * invert_perm(char * perm);

double dist(const Rmatrix & U, const Rmatrix & V);
double dist(const Unitary & U, const Unitary & V);
void init(int n, int m);

struct cmp_hash {
  bool operator()(const hash_t & a, const hash_t & b) const;
};
struct eq_hash {
  bool operator()(const hash_t & a, const hash_t & b) const;
};
struct hasher {
  unsigned int operator()(const hash_t & a) const;
};
hash_t Hash_Unitary(const Unitary & U);
hash_t Hash_Rmatrix(const Rmatrix & U);

void permute(const Rmatrix & U, Rmatrix & V, int i);
void permute_inv(const Rmatrix & U, Rmatrix & V, int i);
void permute_hash(const Rmatrix & U, Rmatrix & V, int i);
void permute_adj_hash(const Rmatrix & U, Rmatrix & V, int i);
bool equiv(const Rmatrix & M, const Rmatrix & N);

Canon canonicalize(const Rmatrix & U, bool sym);

#endif
