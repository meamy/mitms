#include "gate.h"

#include <iostream>
#include <string>
#include <assert.h>
#include <map>
#include <vector>
#include <list>

#define MAX_SEQ 50
#define CLIFF 50
#define SYMMS false

typedef multimap<hash_t, Circuit *, cmp_hash> map_t;
typedef pair    <hash_t, Circuit *>           map_elt;
typedef map_t::iterator                       map_iter;

typedef list<Circuit *>        circuit_list;
typedef circuit_list::iterator circuit_iter;

typedef pair<bool, Circuit *> result;

/* List of base circuits */
circuit_list base_list;

/* Store all unique unitaries up to MAX_SEQ */
map_t circuit_table[MAX_SEQ];
map_t * cliff_temp;

unsigned int numcorrect = 0, numcollision = 0;


/* It seems like this function could be expanded todeal with things like phase, 
   but DO NOT do so -- this way is more efficient */
result find_unitary(hash_t key, Unitary &U, map_t &map) {
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

result find_unitary(hash_t key, Rmatrix &U, map_t &map) {
  pair<map_iter, map_iter> ret;
  map_iter it;
  Rmatrix V(dim, dim);

  ret = map.equal_range(key);
  for (it = ret.first; it != ret.second; ++it) {
    numcollision++;
    (it->second)->to_Rmatrix(V);
    if (U.phase_eq(V)) {
      numcorrect++;
      return result(true, it->second);
    }
  }

  return result(false, NULL);
}

bool generate_cliff(int i) {
  int j, k;
  hash_t key, tmp_key;
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

      flg = false;
      j = 0;
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
        while (!flg && j <= 1) {
          if (j == 0) {
            key = tmp_key = Hash_Rmatrix(V);
          } else {
            V *= phase;
            tmp_key = Hash_Rmatrix(V);
          }

          k = 0;
          while (!flg && k <= i) {
            flg = find_unitary(key, V, cliff_temp[k]).first;
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

/* Generate the iteration set */
void generate_base_circuits(bool cliffords) {
  Circuit * ins;
  int j;

  cout << "Generating iteration set V(n, G)\n";
  if (!cliffords) {
    Gate G;

    // Reset gate G
    for (j = 0; j < num_qubits; j++) {
      G[j] = I;
    }
    ins = new Circuit;
    ins->G = G;
    base_list.push_back(ins);

    /* Generate all the sequences of length i */
    while(!((++G).eye())) {
      ins = new Circuit;
      ins->G = G;
      base_list.push_back(ins);
    }
  } else {
    bool flg = true;
    int i = 0, k, l;
    map_iter it;

    cliff_temp = new map_t[CLIFF];
    while(flg && i < CLIFF) {
      cout << "Cliffords of length " << i + 1 << "\n";
      flg = generate_cliff(i);
      cout << cliff_temp[i].size() << "\n";
      i++;
    }

    if (i >= CLIFF) cout << "ERROR: unique cliffords of length > " << CLIFF << "\n";
    for (j = i-1; j >= 0; j--) {
      for (it = cliff_temp[j].begin(); it != cliff_temp[j].end(); it++) {
        for (k = 0; k < pow(num_qubits, 3); k++) {
          ins = new Circuit;
          for (l = 0; l < num_qubits; l++) {
            if ((k / (int)pow(3, l)) % 3 == 0) {
              ins->G[l] = I;
            } else if ((k / (int)pow(3, l)) % 3 == 1) {
              ins->G[l] = T;
            } else {
              ins->G[l] = Td;
            }
          }
          ins->next = it->second;
          base_list.push_back(it->second);
        }
      }
    }
    delete[] cliff_temp;
  }
  cout << "Iteration set size: " << base_list.size() << "\n";
}

int max (int a, int b) {
  if (a > b) return a;
  else return b;
}

void generate_sequences(int i, circuit_list &L) {
  int j, k;
  bool flg;
  hash_t tmp_key;
  Circuit * tmp_circ, * ins_circ, * tmp;
  Rmatrix V(dim, dim);
  map_iter it;
  circuit_iter c;
  Canon canon_form;
  struct triple * trip;
  Gate G;

  time_t start, end;

  cout << "--------------------------------------\n";
  cout << "Generating sequences of length " << i+1 << "\n" << flush;
  time(&start);

  /* Generate all the sequences of length i */
  for (c = L.begin(); c != L.end(); c++) {
    if (i == 0) {

      /* Create a circuit for the gate */
      tmp_circ = (*c)->append(NULL);

      /* Compute the matrix for circuit C */
      tmp_circ->to_Rmatrix(V);
      canon_form = canonicalize(V, SYMMS);
      // Check to see if it's already synthesized

      while(!canon_form.empty()) {
        trip = &(canon_form.front());
        flg = find_unitary(trip->key, trip->mat, circuit_table[0]).first;
        if (!flg) {
          ins_circ = tmp_circ->permute(trip->permutation);
          if (trip->adjoint) {
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
      G = (*c)->last();
      /* For each circuit of length i ending in gate G */
      for (it = circuit_table[i-1].begin(); it != circuit_table[i-1].end(); it++) {
        flg = false;
        for (j = 0; j < num_qubits; j++) {
          if (!IS_C((it->second)->G[j])) {
            if (G[j] == X && (it->second)->G[j] == X) {
              flg == true;
              for (k = 0; k < num_qubits; k++) {
               flg = flg && !(G[k] == C(j) xor (it->second)->G[k] == C(j));
              }
            } else { 
              flg = flg || (G[j] == adjoint[(it->second)->G[j]]);
            }
          }
        }
        if (!flg) {
          tmp_circ = (*c)->append(it->second);

          tmp_circ->to_Rmatrix(V);
          canon_form = canonicalize(V, SYMMS);

          // Check to see if it's already synthesized
          while(!canon_form.empty()) {
            trip = &(canon_form.front());
            flg = false;
            // Check to see if it's already synthesized in all sequences of length < i
            for (j = 0; j <= i; j++) {
              flg = flg || find_unitary(trip->key, trip->mat, circuit_table[j]).first;
            }
            if (!flg) {
              ins_circ = tmp_circ->permute(trip->permutation);
              if (trip->adjoint) {
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
  time(&end);
  cout << "Time: " << difftime(end, start) << "s\n";
  cout << "# new unitaries: " << circuit_table[i].size() << "\n";
  cout << "--------------------------------------\n" << flush;
}

void exact_search(Rmatrix & U, circuit_list &L) {
  int i, j, k;
  Circuit * ans, * tmp_circ;
  Rmatrix V(dim, dim);
  map_iter it;
  Canon canon_form;
  struct triple * trip;
  int s;

  if (SYMMS) s = 2*num_swaps;
  else s = 1;

  time_t start, end;

  for (i = 0; i < MAX_SEQ; i++) {
    generate_sequences(i, L);
    cout << "Looking for circuits...\n";
    time(&start);
    /* Meet in the middle - Sequences of length 2i + {0, 1} */
    for (j = max(i-1, 0); j <= i; j++) {
      for (it = circuit_table[j].begin(); it != circuit_table[j].end(); it++) {
        for (k = 0; k < s; k++) {

          if (k == 0) {
            tmp_circ = it->second;
          } else if (k % 2 == 0) {
            /* Adjoint case */
            tmp_circ = (it->second)->permute(k/2);
          } else {
            /* Non-adjoint case */
            tmp_circ = ((it->second)->permute(k/2))->adj(NULL);
          }

          tmp_circ->to_Rmatrix(V);
          V *= U;
          canon_form = canonicalize(V, SYMMS);

          while(!canon_form.empty()) {
            trip = &(canon_form.front());
            ans = find_unitary(trip->key, trip->mat, circuit_table[i]).second;
            if (ans != NULL) {
              (it->second)->print(ans);
              cout << "\n" << flush;
            }
            canon_form.pop_front();
          }
          if (k != 0) delete_circuit(tmp_circ);
        }
      }
    }
    time(&end);
    cout << "Time: " << difftime(end, start) << "s\n";
    cout << "equivalent unitary vs equivalent key: " << numcorrect << " / " << numcollision << "\n";
    cout << "--------------------------------------\n" << flush;
  }
}

/*
void approx_search(Unitary & U) {
  hash_t tmp_key;
  int i, j, flg;
  Circuit * ans, * tmp_circ, * tmp;
  Gate G;
  Rmatrix J(dim, dim);
  Unitary V(dim, dim), tmp_mat(dim, dim);
  map_iter it;

  for (i = 0; i < MAX_SEQ; i++) {
    generate_sequences(i);
    // Meet in the middle - Sequences of length 2i + {0, 1} 
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
*/

void find_first_nlr(circuit_list &L) {
  int i, j, k;
  Circuit * ans, * tmp_circ;
  Rmatrix V(dim, dim);
  map_iter it, ti;
  Canon canon_form;
  struct triple * trip;
  int s;

  for (i = 0; i < MAX_SEQ; i++) {
    generate_sequences(i, L);
    cout << "Looking for circuits...\n";
    /* Meet in the middle - Sequences of length 2i + {0, 1} */
    for (j = max(i-1, 0); j <= i; j++) {
      for (it = circuit_table[j].begin(); it != circuit_table[j].end(); it++) {
        for (ti = circuit_table[i].begin(); ti != circuit_table[i].end(); ti++) {
          tmp_circ = (it->second)->append(ti->second);
          tmp_circ->to_Rmatrix(V);
          if (V.is_nonlinear_reversible()) {
            tmp_circ->print();
          }
        }
      }
    }
    cout << "--------------------------------------\n" << flush;
  }
}

int main() {
  init(1, 1);
  map_iter it;
  generate_base_circuits(false);
  numcorrect = 0;
  numcollision = 0;
/*
  Circuit * x = new Circuit;
  x->G[0] = C(2);
  x->G[1] = I; 
  x->G[2] = X; 
  Circuit * y = x->next = new Circuit;
  y->G[0] = I;
  y->G[1] = C(2); 
  y->G[2] = X; 
  Circuit * z = y->next = new Circuit;
  z->G[0] = C(2);
  z->G[1] = C(2); 
  z->G[2] = X; 
*/
  Circuit *h1 = new Circuit(H, new Circuit(T, new Circuit(T, new Circuit(T, NULL))));
  Circuit *h2 = new Circuit(H, new Circuit(Td, new Circuit(H, new Circuit(Td, new Circuit(H, new Circuit(T, h1))))));
  Circuit *h3 = new Circuit(H, new Circuit(T, new Circuit(T, new Circuit(T, h2))));
  Circuit *h4 = new Circuit(H, new Circuit(Td, new Circuit(Td, new Circuit(Td, h3))));
  Circuit *h5 = new Circuit(H, new Circuit(T, new Circuit(T, new Circuit(T, h4))));
  Circuit *h6 = new Circuit(H, new Circuit(T, h5));
  Circuit *h7 = new Circuit(H, new Circuit(T, new Circuit(T, new Circuit(T, h6))));
  Circuit *h8 = new Circuit(H, new Circuit(T, new Circuit(T, new Circuit(T, h7))));
  Circuit *h9 = new Circuit(H, new Circuit(T, new Circuit(T, new Circuit(T, h8))));
  Circuit *fin = new Circuit(X, new Circuit(T, h9));
  Circuit *tst = new Circuit(S, new Circuit(H, new Circuit(Sd, new Circuit(H, NULL))));
  Rmatrix U(dim, dim);
  Unitary V(dim, dim), W(dim, dim);
  W(0, 0) = LaComplex(1, 0);
  W(1, 0) = LaComplex(0, 0);
  W(0, 1) = LaComplex(0, 0);
  W(1, 1) = LaComplex(cos(PI/8), sin(PI/8));
  fin->to_Rmatrix(U);
  U.to_Unitary(V);
  cout << dist(V, W) << "\n";
 // exact_search(U, base_list);

  tst->to_Rmatrix(U);
  cout << V;
//  assert(U.is_nonlinear_reversible());
//  find_first_nlr(base_list);

  return 0;
}
