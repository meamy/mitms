#ifndef CIRCUIT
#define CIRCUIT

#include "gate.h"

/* --------------- Circuits */
/* Circuits are stored in operator order -- i.e. a circuit --X--T-- is stored as T*X */
class Circuit {
  private:
    int steane_cost_helper(int q) const;

  public:
    short unsigned int depth;
    char * circuit;

    Circuit();
    Circuit(int d);
    Circuit(const Circuit & C);

    Circuit & operator=(const Circuit & C);
    gate operator[](int i);
    Circuit & operator++();
    Circuit & cliffpp();

    gate    first() const;
    Circuit next()  const;
    gate    last()  const;

    bool empty() const;
    bool eye() const;

    Circuit copy() const;
    Circuit append(const Circuit & C) const;
    Circuit reverse() const;
    Circuit transform(const char * perm, bool adj) const;
    Circuit transform(int i, bool adj) const;

    void print() const;

    void to_Rmatrix(Rmatrix & U, bool adj) const;
    void to_Rmatrix(Rmatrix & U) const { to_Rmatrix(U, false); }
    void to_Unitary(Unitary & U, bool adj) const;
    void to_Unitary(Unitary & U) const { to_Unitary(U, false); }
    int cost() const;

    void output(ofstream & out, int depth) const;
    void input (ifstream & in, int depth);

    Circuit adj() const 
      { return this->transform(0, true); }
    Circuit permute(const char * perm) const
      { return this->transform(perm, false); }
    Circuit permute(int i) const
      { return this->transform(i, false); }
    Circuit permute_adj(const char * perm) const
      { return this->transform(perm, true); }
    Circuit permute_adj(int i) const
      { return this->transform(i, true); }

    friend Circuit read_circuit(istream & in);
};

void delete_circuit(Circuit & C);
Circuit read_circuit(istream & in);

void test_circuit();

#endif
