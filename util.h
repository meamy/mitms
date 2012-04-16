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
#define C(x)          (x | 0x80)
#define IS_C(x)       (x & 0x80)
#define GET_TARGET(x) (x & 0x7F)

/* Definitions */
#define SUBSPACE_SIZE 1   // Size of subspace we store
#define PRECISION 1       // Precision
#define PHASE             // Whether we mod out phase
#define MAX_SEQ 50        // Max circuit length
#define CLIFF 50          // Max depth we try to find clifford circuits
#define SYMMS true        // Whether we mod out symmetries
#define CHECK_EQUIV false // Whether we check to make sure two circuits are equiv



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
typedef LaGenMatComplex hash_t;

struct triple {
  Rmatrix mat;
  hash_t key;
  bool adjoint;
  int permutation;
};

int fac(int n);
int to_lexi(char * perm);
char * from_lexi(int n);

double dist(const Rmatrix & U, const Rmatrix & V);
double dist(const Unitary & U, const Unitary & V);
void init(int n, int m);

struct cmp_hash {
  bool operator()(const hash_t & a, const hash_t & b);
};
bool operator<(const hash_t & a, const hash_t & b);
bool operator==(const hash_t & a, const hash_t & b);
hash_t Hash_Unitary(const Unitary & U);
hash_t Hash_Rmatrix(const Rmatrix & U);

unsigned int hasher(hash_t & U);

void permute(const Rmatrix & U, Rmatrix & V, int i);

Canon canonicalize(const Rmatrix & U, bool sym);

#endif
