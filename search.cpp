#include "gate.h"

#include <iostream>
#include <string>
#include <assert.h>
#include <map>
#include <vector>

#define MAX_SEQ 50
#define CLIFF 50
//define SYMMETRIES

typedef multimap<hash_t, Circuit *, cmp_hash> map_t;
typedef     pair<hash_t, Circuit *> map_elt;
typedef map_t::iterator          map_iter;

typedef     pair<bool,   Circuit *> result;

map_t * cliff_temp;
map_t cliffords;

/* Store all unique unitaries up to MAX_SEQ */
map_t circuit_table[MAX_SEQ];


/* It seems like this function could be expanded todeal with things like phase, 
   but DO NOT do so -- this way is more efficient */
result find_unitary(hash_t key, Unitary &U, map_t map) {
  pair<map_iter, map_iter> ret;
  map_iter it, ti;
  Unitary V(dim, dim);

  ret = map.equal_range(key);
  for (it = ret.first; it != ret.second; ++it) {
    (it->second)->to_Unitary(V);
    if (dist(U, V) < 0.01) {
      return result(true, it->second);
    }
  }

  return result(false, NULL);
}

result find_unitary(hash_t key, Rmatrix &U, map_t map) {
  pair<map_iter, map_iter> ret;
  map_iter it;
  Rmatrix V(dim, dim);

  ret = map.equal_range(key);
  for (it = ret.first; it != ret.second; ++it) {
    (it->second)->to_Rmatrix(V);
    if (U.phase_eq(V)) {
      return result(true, it->second);
    }
  }

  return result(false, NULL);
}

bool generate_cliff(int i) {
  int j, k;
  hash_t tmp_key, key;
  Gate G;
  Circuit * tmp_circ;
  Rmatrix V(dim, dim);
  map_iter it;
  bool ret = false, flg = false;
  Elt phase(0, 1, 0, 0, 0);

  // Reset gate G
  for (j = 0; j < num_qubits; j++) {
    G[j] = I;
  }
  if (i == 0) {
    tmp_circ = new Circuit;
    tmp_circ->G = G;
    tmp_circ->to_Rmatrix(V);
    cliff_temp[0].insert(map_elt(Hash_Rmatrix(V), tmp_circ));
  }

  /* Generate all the sequences of length i */
  while(!((G.cliffpp()).eye())) {
    if (i == 0) {

      /* Create a circuit for the gate */
      tmp_circ = new Circuit;
      tmp_circ->G = G;
      tmp_circ->next = NULL;
      tmp_circ->to_Rmatrix(V);

      /* Compute the matrix for circuit C */
      j = 0;
      flg = false;
      while(!flg && j < 8) {
        if (j == 0) {
          key = tmp_key = Hash_Rmatrix(V);
        } else {
          V *= phase;
          tmp_key = Hash_Rmatrix(V);
        }

        flg = find_unitary(tmp_key, V, cliff_temp[i]).first;
        j++;
      }
      if (!flg) {
        cliff_temp[i].insert(map_elt(key, tmp_circ));
        ret = true;
      } else {
        delete tmp_circ;
      }
    } else {
      /* For each circuit of length i ending in gate G */
      for (it = cliff_temp[i-1].begin(); it != cliff_temp[i-1].end(); it++) {
        tmp_circ = new Circuit;
        tmp_circ->G = G;
        tmp_circ->next = it->second;

        tmp_circ->to_Rmatrix(V);

        j = 0;
        flg = false;
        /* For each phase, try to find it. If we don't find it for any phase, insert */
        while(!flg && j < 8) {
          if (j == 0) {
            key = tmp_key = Hash_Rmatrix(V);
          } else {
            V *= phase;
            tmp_key = Hash_Rmatrix(V);
          }

          k = 0;
          while(!flg && k <= i) {
            flg = find_unitary(tmp_key, V, cliff_temp[k]).first;
            k++;
          }
          j++;
        }

        if (!flg) {
          cliff_temp[i].insert(map_elt(key, tmp_circ));
          ret = true;
        } else {
          delete tmp_circ;
        }
      }
    }
  }
  return ret;
}

void generate_cliffords() {
  bool flg = true;
  int i = 0, j;
  map_iter it;

  cliff_temp = new map_t[CLIFF];
  while(flg && i < CLIFF) {
    cout << "Cliffords of length " << i + 1 << "\n";
    flg = generate_cliff(i++);
  }
  if (i >= CLIFF) cout << "ERROR: unique cliffords of length > " << CLIFF << "\n";
  for (j = 0; j < i; j++) {
    for (it = cliff_temp[j].begin(); it != cliff_temp[j].end(); it++) {
      cliffords.insert(*it);
    }
  }
  delete[] cliff_temp;
}

int max (int a, int b) {
  if (a > b) return a;
  else return b;
}

#ifndef SYMMETRIES
void generate_sequences(int i) {
  bool flg;
  int j, k;
  hash_t tmp_key, key;
  Gate G;
  Circuit * tmp_circ;
  Rmatrix V(dim, dim);
  map_iter it;
  Elt phase(0, 1, 0, 0, 0);

  cout << "Generating sequences of length " << i+1 << "\n" << flush;
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

      j = 0;
      flg = false;
      while(!flg && j < 8) {
        if (j == 0) {
          key = tmp_key = Hash_Rmatrix(V);
        } else {
          V *= phase;
          tmp_key = Hash_Rmatrix(V);
        }

        flg = find_unitary(tmp_key, V, circuit_table[i]).first;
        j++;
      }
      if (!flg) {
        circuit_table[i].insert(map_elt(key, tmp_circ));
      } else {
        delete tmp_circ;
      }

    } else {
      /* For each circuit of length i ending in gate G */
      for (it = circuit_table[i-1].begin(); it != circuit_table[i-1].end(); it++) {
        tmp_circ = new Circuit;
        tmp_circ->G = G;
        tmp_circ->next = it->second;

        tmp_circ->to_Rmatrix(V);

        j = 0;
        flg = false;
        /* For each phase, try to find it. If we don't find it for any phase, insert */
        while(!flg && j < 8) {
          if (j == 0) {
            key = tmp_key = Hash_Rmatrix(V);
          } else {
            V *= phase;
            tmp_key = Hash_Rmatrix(V);
          }

          k = 0;
          while(!flg && k <= i) {
            flg = find_unitary(tmp_key, V, circuit_table[k]).first;
            k++;
          }
          j++;
        }

        if (!flg) {
          circuit_table[i].insert(map_elt(key, tmp_circ));
        } else {
          delete tmp_circ;
        }
      }
    }
  }
}

void exact_search(Rmatrix & U) {
  hash_t tmp_key;
  int i, j, flg;
  Circuit * ans, * tmp_circ, * tmp;
  Gate G;
  Rmatrix V(dim, dim), tmp_mat(dim, dim);
  map_iter it;

  /* Insert the identity */
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
    generate_sequences(i);
    /* Meet in the middle - Sequences of length 2i + {0, 1} */
    for (j = max(i-1, 0); j <= i; j++) {
      for (it = circuit_table[j].begin(); it != circuit_table[j].end(); it++) {
        tmp_circ = it->second;
        tmp_circ->to_Rmatrix(tmp_mat);
        tmp_mat.adj(V);
        V*=U;
        tmp_key = Hash_Rmatrix(V);
        ans = find_unitary(tmp_key, V, circuit_table[i]).second;
        if (ans != NULL) {
          (it->second)->print(ans);
          cout << "\n" << flush;
        }
      }
    }
  }
}

#else
void generate_sequences(int i) {
  int j, flg;
  hash_t tmp_key;
  Gate G;
  Circuit * tmp_circ, * ins_circ, * tmp;
  Rmatrix V(dim, dim);
  map_iter it;
  Canon canon_form;
  struct triple * trip;

  cout << "Generating sequences of length " << i+1 << "\n" << flush;
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
      // Check to see if it's already synthesized

      while(!canon_form.empty()) {
        trip = &(canon_form.front());
        flg = find_unitary(trip->key, trip->mat, circuit_table[0]).first;
        if (flg == 0) {
          ins_circ = tmp_circ->permute(trip->permutation);
          if (trip->adjoint != 0) {
            tmp = ins_circ;
            ins_circ = tmp->adj(NULL);
            delete_circuit(tmp);
          }
          circuit_table[0].insert(map_elt(trip->key, ins_circ));
        }
        canon_form.pop_front();
      }
      delete tmp_circ;
    } else {
      /* For each circuit of length i ending in gate G */
      for (it = circuit_table[i-1].begin(); it != circuit_table[i-1].end(); it++) {
        tmp_circ = new Circuit;
        tmp_circ->G = G;
        tmp_circ->next = it->second;

        tmp_circ->to_Rmatrix(V);
        canon_form = canonicalize(V);

        // Check to see if it's already synthesized
        while(!canon_form.empty()) {
          trip = &(canon_form.front());
          flg = 0;
          // Check to see if it's already synthesized in all sequences of length < i
          for (j = 0; j <= i; j++) {
            flg = flg || find_unitary(trip->key, trip->mat, circuit_table[j]).first;
          }
          if (flg == 0) {
            ins_circ = tmp_circ->permute(trip->permutation);
            if (trip->adjoint != 0) {
              tmp = ins_circ;
              ins_circ = tmp->adj(NULL);
              delete_circuit(tmp);
            }
            circuit_table[i].insert(map_elt(trip->key, ins_circ));
          }
          canon_form.pop_front();
        }
        delete tmp_circ;
      }
    }
  }
}

void exact_search(Rmatrix & U) {
  hash_t tmp_key;
  int i, j, k, flg;
  Circuit * ans, * tmp_circ, * tmp;
  Gate G;
  Rmatrix V(dim, dim), tmp_mat(dim, dim);
  Rmatrix syms[num_swaps][2];
  map_iter it;

  /* Find all symmetries */
  for (i = 0; i < num_swaps; i++) {
    permute(U, syms[i][0], i);
    syms[i][0].adj(syms[i][1]);
  }

  /* Insert the identity */
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
    generate_sequences(i);
    /* Meet in the middle - Sequences of length 2i + {0, 1} */
    for (j = max(i-1, 0); j <= i; j++) {
      for (it = circuit_table[j].begin(); it != circuit_table[j].end(); it++) {
        for (k = 0; k < 2*num_swaps; k++) {
          tmp_circ = it->second;
          tmp_circ->to_Rmatrix(tmp_mat);
          tmp_mat.adj(V);
          V*=syms[k/2][k%2];
          tmp_key = Hash_Rmatrix(V);
          ans = find_unitary(tmp_key, V, circuit_table[i]).second;
          if (ans != NULL) {
            (it->second)->print(ans);
            cout << "\n" << flush;
          }
        }
      }
    }
  }
}
#endif

void approx_search(Unitary & U) {
  hash_t tmp_key;
  int i, j, flg;
  Circuit * ans, * tmp_circ, * tmp;
  Gate G;
  Rmatrix J(dim, dim);
  Unitary V(dim, dim), tmp_mat(dim, dim);
  map_iter it;

  /* Insert the identity */
  for (j = 0; j < num_qubits; j++) {
    G[j] = I;
  }
  tmp_circ = new Circuit;
  tmp_circ->G = G;
  tmp_circ->next = NULL;
  tmp_circ->to_Rmatrix(J);
  tmp_key = Hash_Rmatrix(J);
  circuit_table[0].insert(map_elt(tmp_key, tmp_circ));


  for (i = 0; i < MAX_SEQ; i++) {
    generate_sequences(i);
    /* Meet in the middle - Sequences of length 2i + {0, 1} */
    for (j = max(i-1, 0); j <= i; j++) {
      for (it = circuit_table[j].begin(); it != circuit_table[j].end(); it++) {
        tmp_circ = it->second;
        tmp_circ->to_Unitary(tmp_mat);
        Blas_Mat_Mat_Mult(tmp_mat, U, V, true, false);
        tmp_key = Hash_Unitary(V);
        ans = find_unitary(tmp_key, V, circuit_table[i]).second;
        if (ans != NULL) {
          (it->second)->print(ans);
          cout << "\n" << flush;
        }
      }
    }
  }
}

int main() {
  init(3);
  map_iter it;
  /*
  generate_cliffords();
  for (it = cliffords.begin(); it != cliffords.end(); it++) {
    (it->second)->print();
    cout << "\n";
  }
  */
  Circuit * x = new Circuit;
  x->G[0] = C(2);
  x->G[1] = C(2); 
  x->G[2] = X; 
  /*
  Unitary U = Unitary::zeros(dim, dim);
  U(0, 0) = LaComplex(1, 0);
  U(1, 1) = LaComplex(cos(M_PI/8), sin(M_PI/8));

  x->G[2] = X; 
  Circuit * y = x->next = new Circuit;
  y->G[0] = C(1);
  y->G[1] = X; 
  y->G[2] = C(1); 
  Circuit * z = y->next = new Circuit;
  z->G[0] = C(2);
  z->G[1] = C(2); 
  z->G[2] = X; 
  */
  Rmatrix U(dim, dim);
  x->to_Rmatrix(U);
  exact_search(U);
  /*
  approx_search(U);
  */

  return 0;
}
