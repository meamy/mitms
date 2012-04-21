#ifndef GATE
#define GATE

#include "util.h"

/* -------------- Gates */
class Gate {
  private:
    char * gates;

    void tensor(Rmatrix & U) const;
    void tensor(Rmatrix & U, bool adj) const;
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
    void to_Rmatrix(Rmatrix & U, bool adj) const;
    void to_Unitary(Unitary & U) const;
    void permute(Gate & G, char *  perm) const;
    void permute(Gate & G, int i) const;
    void permute_adj(Gate & G, char * perm) const;
    void permute_adj(Gate & G, int i) const;
    void print() const;
};

struct gate_eq {
  bool operator()(const Gate & A, const Gate & B) const { return A == B; }
};

struct gate_hasher {
  unsigned int operator()(const Gate & R) const;
};

void init_ht();

void gate_test();

#endif
