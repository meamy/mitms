#ifndef UTIL
#define UTIL

#include "circuit.h"
#include <list>
#include <iostream>
#include <string>
#include <assert.h>
#include <map>
#include <limits>
#include <fstream>

#include <gmc.h>
#include <blas3pp.h>
#include <blas2pp.h>
#include <blaspp.h>
#include <laslv.h>
#include <lavc.h>

/* Configs */
#define SUBSPACE_SIZE 1   // Size of subspace we store
#define SUBSPACE_ABS  false// take the absolute value of the subspace matrix
#define PRECISION 1       // Precision
#define PHASE false       // Whether we mod out phase
#define MAX_SEQ 50        // Max circuit length
#define CLIFF 50          // Max depth we try to find clifford circuits
#define SYMMS true        // Whether we mod out symmetries
#define CHECK_EQUIV true  // Whether we check to make sure two circuits are equiv
#define ORDERED true      // Whether we should use an ordered map
#define TENSORS true      // Whether to store gates as tensor products of gates
#define TDEPTH  true      // Whether we want to search by T-depth
#define SERIALIZE false    // Whether to serialize the generated circuits


using namespace std;

typedef LaGenMatComplex Unitary;
typedef list< struct triple > Canon;
 
typedef LaGenMatComplex hash_t;
typedef LaGenMatComplex subs_t;

struct triple {
  Rmatrix mat;
  hash_t key;
  bool adjoint;
  int permutation;

  triple() {}
  triple(Rmatrix m, hash_t k, bool adj, int perm) { mat = m; key = k; adjoint = adj; permutation = perm; }
};

int max (int a, int b);

double dist(const Rmatrix & U, const Rmatrix & V);
double dist(const Unitary & U, const Unitary & V);
void init_util();

bool operator==(const hash_t & a, const hash_t & b);

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

Canon * canonicalize(const Rmatrix & U, bool sym);
Canon * canonicalize(const Rmatrix & U);

void output_key(ofstream & out, const hash_t & key);
void input_key (ifstream & in, hash_t & key);

#endif
