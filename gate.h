#ifndef GATE
#define GATE

#include "matrix.h"

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
           void to_Rmatrix(Rmatrix & U, bool adj) const;
    inline void Gate::to_Rmatrix(Rmatrix & U) const { this->to_Rmatrix(U, false); }
    void to_Unitary(Unitary & U) const;
    void permute(Gate & G, const char *  perm) const;
    void permute(Gate & G, int i) const;
    void permute_adj(Gate & G, const char * perm) const;
    void permute_adj(Gate & G, int i) const;
    void print() const;

    void output(ofstream & out) const;
    void input(ifstream & in);
};

struct gate_eq {
  bool operator()(const Gate & A, const Gate & B) const { return A == B; }
};

bool nontrivial_id(const Gate & A, const Gate & B);

unsigned int gate_hasher(const Gate & R);

void init_gate();

void test_gate();

#endif
