#ifndef GATE
#define GATE

#define LA_COMPLEX_SUPPORT

#include <iostream>
#include <string>
#include <assert.h>
#include <map>
#include <list>
#include <limits>

#include <gmc.h>
#include <blas3pp.h>
#include <blas2pp.h>
#include <blaspp.h>
#include <laslv.h>
#include <lavc.h>

using namespace std;

typedef LaGenMatComplex Unitary;

extern int num_qubits;

/* -------------- Gates */
class Gate {
  private:
    char * gates;

    void tensor(Unitary & U) const;
  public:
    Gate();
    Gate(const Gate & G);
    ~Gate();

    char & operator[](int i) const;
    Gate & operator=(const Gate & G);

    void adj(Gate & G) const;
    void to_Unitary(Unitary & U) const;
    void permute(Gate & G, char *  perm) const;
    void print() const;
};

/* --------------- Circuits */
class Circuit {
  public:
    Gate G;
    Circuit * next;

    void print_circuit() const;
    Circuit * adj(Circuit * last) const;
    void to_Unitary(Unitary & U) const;
    void print() const;
};

/*
typedef struct Circuit {
	Gate G;
	struct Circuit * next;
} Circuit;
*/
typedef list< pair<char, int> > Canon;

void print_circuit(const Circuit * C);
Unitary Unitary_of_Circuit(const Circuit * C);

double dist(const Unitary & U, const Unitary & V);
void init(int n);

int Hash_Unitary(const Unitary &U);
Canon canonicalize(const Unitary & U);
void test();

#endif
