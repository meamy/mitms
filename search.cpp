#include "gate.h"

#include <iostream>
#include <string>
#include <assert.h>
#include <map>

#define MAX_SEQ 16

typedef multimap<double, Circuit *> map_t;
typedef     pair<double, Circuit *> map_elt;
typedef map_t::iterator          map_iter;


/* Store all unique unitaries up to MAX_SEQ */
map_t circuit_table[MAX_SEQ];

map_elt find_unitary(int key, Unitary &U, int l) {
  pair<map_iter, map_iter> ret;
  map_iter it, ti;
  Unitary V(dim, dim);

  ret = circuit_table[l].equal_range((double)key);
  for (it = ret.first; it != ret.second; ++it) {
    (it->second)->to_Unitary(V);
    if (dist(U, V) < 0.0001) {
      return map_elt(1, it->second);
    }
  }

  return map_elt(0, NULL);
}

map_elt find_unitary(double key, Rmatrix &U, int l) {
  pair<map_iter, map_iter> ret;
  map_iter it;
  Rmatrix V(dim, dim);

  ret = circuit_table[l].equal_range(key);
  for (it = ret.first; it != ret.second; ++it) {
    (it->second)->to_Rmatrix(V);
    if (U.phase_eq(V)) {
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


void exact_search(Rmatrix & U) {
  double key = Hash_Rmatrix(U), tmp_key;
  int i, j, flg;
  Circuit * ans, * tmp_circ;
  Gate G;
  Rmatrix V(dim, dim), tmp_mat(dim, dim);
  map_iter it;
  Canon canon_form;

  for (j = 0; j < num_qubits; j++) {
    G[j] = I;
  }
  tmp_circ = new Circuit;
  tmp_circ->G = G;
  tmp_circ->next = NULL;
  tmp_circ->to_Rmatrix(V);
  tmp_key = Hash_Rmatrix(V);
  circuit_table[0].insert(map_elt(tmp_key, tmp_circ));

  for (i = 0; i < MAX_SEQ; i++) {
    cout << "Generating sequences of length " << i+1 << "\n";
    // Reset gate G
    for (j = 0; j < num_qubits; j++) {
      G[j] = I;
    }

    /* Generate all the sequences of length i */
    while(!((++G).eye())) {
      if (i == 0) {
        /* Create a circuit for the gate */
        tmp_circ = new Circuit;
        tmp_circ->G = G;
        tmp_circ->next = NULL;

        /* Compute the matrix for circuit C */
        tmp_circ->to_Rmatrix(V);
        canon_form = canonicalize(V);
        /*
        while(!canon_form.empty()) {
          trip = &(canon_form.front());
          flg = find_unitary(trip->key, trip->mat, 0).first;
          if (flg == 0) {
            ins_circ = tmp_circ->permute(trip->permutation);
            if (trip->adjoint != 0) {
              tmp = ins_circ;
              ins_circ = tmp->adj();
              tmp->full_delete();
            }
            circuit_table[0].insert(map_elt(trip->key, ins_circ));
          }
        }
        delete tmp_circ;
*/

        tmp_key = Hash_Rmatrix(V);

        /* Check to see if it's already synthesized */
        flg = find_unitary(tmp_key, V, 0).first;
        if (flg == 0) circuit_table[0].insert(map_elt(tmp_key, tmp_circ));
        else delete tmp_circ;
      } else {
        /* For each circuit of length i ending in gate G */
        for (it = circuit_table[i-1].begin(); it != circuit_table[i-1].end(); it++) {
          tmp_circ = new Circuit;
          tmp_circ->G = G;
          tmp_circ->next = it->second;

          tmp_circ->to_Rmatrix(V);
          canon_form = canonicalize(V);
          tmp_key = Hash_Rmatrix(V);
          flg = 0;
          /* Check to see if it's already synthesized in all sequences of length < i */
          for (j = 0; j <= i; j++) {
            flg = flg || find_unitary(tmp_key, V, j).first;
          }
          if (flg == 0) circuit_table[i].insert(map_elt(tmp_key, tmp_circ));
          else delete tmp_circ;
        }
      }
    }
    /* Meet in the middle */
    /* Sequences of length 2i - 1 */
    if (i > 0) {
      for (it = circuit_table[i-1].begin(); it != circuit_table[i-1].end(); it++) {
        tmp_circ = it->second;
        tmp_circ->to_Rmatrix(tmp_mat);
        tmp_mat.adj(V);
        V*=U;
        tmp_key = Hash_Rmatrix(V);
        ans = find_unitary(tmp_key, V, i).second;
        if (ans != NULL) {
          (it->second)->print(ans);
          cout << "\n";
        }
      }
    }
    /* Sequences of length 2i */
    for (it = circuit_table[i].begin(); it != circuit_table[i].end(); it++) {
      tmp_circ = it->second;
      tmp_circ->to_Rmatrix(tmp_mat);
      tmp_mat.adj(V);
      V*=U;
      tmp_key = Hash_Rmatrix(V);
      ans = find_unitary(tmp_key, V, i).second;
      if (ans != NULL) {
        (it->second)->print(ans);
        cout << "\n";
      }
    }
  }
}


int main() {
  init(2);
  Circuit * x = new Circuit;
  x->G[0] = C(1);
  x->G[1] = H;
  char tst[] = {1, 0};
  Circuit * y = x->permute(tst);
  x->print();
  y->print();
  Rmatrix U(dim, dim);
  x->to_Rmatrix(U);
  exact_search(U);

  return 0;
}
