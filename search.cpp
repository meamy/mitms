#include "gate.h"

#include <iostream>
#include <string>
#include <assert.h>
#include <map>

#define MAX_SEQ 16

typedef multimap<int, Circuit *> map_t;
typedef     pair<int, Circuit *> map_elt;
typedef map_t::iterator          map_iter;


/* Store all unique unitaries up to length 16 */
map_t circuit_table[MAX_SEQ];

map_elt find_unitary(int key, Unitary &U, int l) {
  pair<map_iter, map_iter> ret;
  map_iter it, ti;
  Unitary V(dim, dim);

  ret = circuit_table[l].equal_range(key);
  for (it = ret.first; it != ret.second; ++it) {
    (it->second)->to_Unitary(V);
    if (dist(U, V) < 0.0001) {
      return map_elt(1, it->second);
    }
  }

  return map_elt(0, NULL);
}

map_elt find_unitary(int key, Rmatrix &U, int l) {
  pair<map_iter, map_iter> ret;
  map_iter it, ti;
  Rmatrix V(dim, dim);

  ret = circuit_table[l].equal_range(key);
  for (it = ret.first; it != ret.second; ++it) {
    (it->second)->to_Rmatrix(V);
    if (V.phase_eq(V)) {
      return map_elt(1, it->second);
    }
  }

  return map_elt(0, NULL);
}

/* For each qubit, for each possible gate at that qubit we recurse */
void choose_gate(Gate &G, int q, int rem, int n);
void choose_gate(Gate &G, int q, int rem, int n) {
  int i;

  if (rem == 0) {
    pair<map_iter, map_iter> ret;
    map_iter it, ti;
    Circuit * C;
    int key;
    char flg = 0;

    if (n == 0) {
      C = new Circuit;
      C->next = NULL;
      for (i = 0; i < num_qubits; i++) {
        C->G[i] = G[i];
      }

      Rmatrix U(dim, dim);
      C->to_Rmatrix(U);

      key = Hash_Rmatrix(U);

      flg = find_unitary(key, U, n).first;
      if (flg != 1) circuit_table[n].insert(map_elt(key, C));
      else delete C;

    } else {
      for (ti = circuit_table[n-1].begin(); ti != circuit_table[n-1].end(); ti++) {
        C = new Circuit;
        C->next = ti->second;
        for (i = 0; i < num_qubits; i++) {
          C->G[i] = G[i];
        }

        Rmatrix U(dim, dim);
        C->to_Rmatrix(U);
        key = Hash_Rmatrix(U);
        flg = 0;
        for (int j = 0; j <= n; j++) {
          flg = flg || find_unitary(key, U, j).first;
        }
        if (flg == 0) circuit_table[n].insert(map_elt(key, C));
        else delete C;
      }
    }
  } else if (!IS_C(G[q])) {
    /* Try each single qubit gate in position q */
    for (i = 0; i < 8; i++) {
      G[q] = i;
      choose_gate(G, q+1, rem-1, n);
    }
    /* Try each possible CNOT in position q */
    if (rem > 1) {
      for (i = q+1; i < num_qubits; i++) {
        if (!IS_C(G[i])) { 
          G[q] = C(i);
          G[i] = X;
          choose_gate(G, q+1, rem-2, n);

          G[q] = X;
          G[i] = C(q);
          choose_gate(G, q+1, rem-2, n);

          G[q] = 0;
          G[i] = 0;
        }
      }
    }
  } else {
    choose_gate(G, q+1, rem, n);
  }
}


void do_it(Unitary &U) {
  int key = Hash_Unitary(U);
  int i, j, tmp_key;
  Circuit * ans, * tmp_circ;
  Gate G;
  Unitary t, tmp_unit(dim, dim);
  map_iter it;

  for (i=0; i<16; i++) {
    cout << "Phase " << i << "\n";
    for (j=0; j<num_qubits; j++) {
      G[j] = 0;
    }
    choose_gate(G, 0, num_qubits, i);

    /* Check to see if we've synthesized it */
    ans = find_unitary(key, U, i).second;
    if (ans != NULL) {
      ans->print();
      return;
    }

    /* Meet in the middle stuff */
    if (i > 0) {
      for (it = circuit_table[i-1].begin(); it != circuit_table[i-1].end(); it++) {
        tmp_circ = (it->second)->adj(NULL);
        tmp_circ->to_Unitary(t);
        Blas_Mat_Mat_Mult(t, U, tmp_unit, false, false, 1, 0);
        tmp_key = Hash_Unitary(tmp_unit);
        ans = find_unitary(tmp_key, tmp_unit, i).second;
        if (ans != NULL) {
          (it->second)->print();
          ans->print();
        }
      }
    }

    for (it = circuit_table[i].begin(); it != circuit_table[i].end(); it++) {
      tmp_circ = (it->second)->adj(NULL);
      tmp_circ->to_Unitary(t);
      Blas_Mat_Mat_Mult(t, U, tmp_unit, false, false, 1, 0);
      tmp_key = Hash_Unitary(tmp_unit);
      ans = find_unitary(tmp_key, tmp_unit, i).second;
      if (ans != NULL) {
        (it->second)->print();
        ans->print();
      }
    }

  }
}


int main() {
  init(3);
  return 0;
}
