#ifndef GATE
#define GATE

#include <iostream>
#include <string>
#include <assert.h>
#include <map>
#include <list>
#include <limits>

#include "matrix.h"
#include <blas3pp.h>
#include <blas2pp.h>
#include <blaspp.h>
#include <laslv.h>
#include <lavc.h>

#define PI 3.14159
/* Number of basis gates */
#define basis_size 9

/* Clifford group on 1 qubit */
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

/* Size of the subspace */
#define SUBSPACE_SIZE 1
#define PRECISION 1


using namespace std;

typedef LaGenMatComplex Unitary;
typedef list< struct triple > Canon;
typedef LaGenMatComplex hash_t;

extern int num_qubits;
extern int num_swaps;
extern int dim;
extern int reduced_dim;
extern const char adjoint[];

/* -------------- Gates */
class Gate {
  private:
    char * gates;

    void tensor(Rmatrix & U) const;
    bool valid_gate();
    void increment(bool t);
  public:
    Gate();
    Gate(const Gate & G);
    ~Gate();

    char & operator[](int i) const;
    Gate & operator=(const Gate & G);
    Gate & operator++();
    Gate & cliffpp();
    const bool operator==(const Gate & G) const;
    const bool operator!=(const Gate & G) const { return !(*this == G); }

    const bool eye() const;
    void adj(Gate & G) const;
    void to_Rmatrix(Rmatrix & U) const;
    void to_Unitary(Unitary & U) const;
    void permute(Gate & G, char *  perm) const;
    void permute(Gate & G, int i) const;
    void permute_adj(Gate & G, char * perm) const;
    void permute_adj(Gate & G, int i) const;
    void print() const;
};

/* --------------- Circuits */
class Circuit {
  public:
    Gate G;
    Circuit * next;
    Circuit();
    Circuit(char g, Circuit * next);

    Circuit * reverse(Circuit * last) const;
    Circuit * adj(Circuit * last) const;
    Circuit * permute(char * perm) const;
    Circuit * permute(int i) const;
    Circuit * permute_adj(char * perm, Circuit * last) const;
    Circuit * permute_adj(int i, Circuit * last) const;
    Circuit * append(Circuit * C) const;
    const Gate & last() const;
    void to_Rmatrix(Rmatrix & U) const;
    void to_Unitary(Unitary & U) const;
    void print() const;
    void print(Circuit * snd) const;
};

void delete_circuit(Circuit * circ);

struct triple {
  Rmatrix mat;
  hash_t key;
  bool adjoint;
  int permutation;
};

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

void permute(const Rmatrix & U, Rmatrix & V, int i);

Canon canonicalize(const Rmatrix & U, bool sym);

void gate_test();
void circuit_test();

#endif
