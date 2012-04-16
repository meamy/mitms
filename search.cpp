#include "gate.h"

#include <iostream>
#include <string>
#include <assert.h>
#include <map>
#include <vector>
#include <list>

#define MAX_SEQ 50
#define CLIFF 50
#define SYMMS true
#define CHECK_EQUIV false

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
/* Store all unique unitaries with key given by their mapping to the |0> subspace */
map_t left_table[MAX_SEQ];
map_t * cliff_temp;

unsigned int numcorrect = 0, numcollision = 0;

/* Determine if there is a nontrivial pair of gates that multiply to the identity */
bool nontrivial_id(const Gate & G, const Circuit * circ) {
  int i, j;
  bool ret = false;

  for (i = 0; i < num_qubits; i++) {
    if (!IS_C(circ->G[i])) {
      if (G[i] == X && circ->G[i] == X) {
        ret = true;
        for (j = 0; j < num_qubits; j++) {
          ret = ret && !(G[j] == C(i) xor circ->G[j] == C(i));
        }
      } else { 
        ret = ret || (G[i] == adjoint[circ->G[i]]);
      }
    }
  }
  return ret;
}

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
  Rmatrix V(dim, dim), W(U.rows(), U.cols()), *tmp;

  /* Find the range of elements with equal keys */
  ret = map.equal_range(key);
  for (it = ret.first; it != ret.second; ++it) {
    numcollision++;
    (it->second)->to_Rmatrix(V);
    /* Determine what submatrix we're comparing */
    if (U.rows() != dim || U.cols() != dim) {
      /* In general, if U has different dimension and we found it in map, the keys
         in map should refer to this submatrix of the actual circuit store */
      V.submatrix(0, 0, U.rows(), U.cols(), W);
      tmp = &W;
    } else {
      tmp = &V;
    }
    /* If we want to play it safe, do this check */
    if (!CHECK_EQUIV || U.phase_eq(*tmp)) {
      numcorrect++;
      return result(true, it->second);
    }
  }

  return result(false, NULL);
}

/* Insert a circuit in a forest of trees, up to a specific depth */
bool insert_tree(Circuit * circ, map_t * tree_list, int depth, bool symm) {
  Rmatrix V(dim, dim), W(reduced_dim, dim);
  Canon canon_form;
  struct triple * trip;
  bool flg = false, ret = false, del = false;
  int i;
  Circuit * ins_circ;

  /* Compute the matrix for circuit C */
  circ->to_Rmatrix(V);
  /* Compute the canonical form */
  canon_form = canonicalize(V, symm);

  // Check to see if it's already synthesized
  // We only need to check the canonical form to get all permutations, inversions and phases
  // Also, we store all canonical forms so that when we search, we only have to look at one canonical form
  while (!canon_form.empty()) {
    trip = &(canon_form.front());
    // Search in each tree up to depth
    flg = false;
    for (i = 1; i <= depth; i++) {
      flg = flg || find_unitary(trip->key, trip->mat, tree_list[i]).first;
    }
    // If none of the searches turned up anything...
    if (!flg) {
      // Hacky if structure, but will be more efficient in the end
      // This would be a great time for a match statement...
      if (trip->permutation == 0) {
        if (trip->adjoint == 0) {
          ins_circ = circ; 
          del = true;
        } else {
          ins_circ = circ->adj(NULL);
        }
      } else {
        if (trip->adjoint == 0) {
          ins_circ = circ->permute(trip->permutation);
        } else {
          ins_circ = circ->permute_adj(trip->permutation, NULL);
        }
      }
      tree_list[depth].insert(map_elt(trip->key, ins_circ));
      ret = true;

      // Searching for reduced dimensions
      if (dim != reduced_dim) {
        (trip->mat).submatrix(0, 0, reduced_dim, dim, W);
        left_table[i].insert(map_elt(Hash_Rmatrix(W), ins_circ));
      }
      canon_form.pop_front();
    } else {
      canon_form.clear();
    }
  }

  if (!del) delete circ;
  return ret;
}


bool generate_cliff(int i) {
  int j, k;
  Gate G;
  Circuit * tmp_circ;
  Rmatrix V(dim, dim);
  map_iter it;
  bool ret = false, flg = false;

  // Reset gate G
  for (j = 0; j < num_qubits; j++) {
    G[j] = I;
  }
  if (i == 0) {
    G.to_Rmatrix(V);
    cliff_temp[i].insert(map_elt(Hash_Rmatrix(V), NULL));
    return true;
  }
  else if (i == 1) {
    tmp_circ = new Circuit;
    tmp_circ->G = G;
    tmp_circ->to_Rmatrix(V);
    cliff_temp[i].insert(map_elt(Hash_Rmatrix(V), tmp_circ));
  }

  /* Generate all the sequences of length i */
  while(!((G.cliffpp()).eye())) {
    /* For each circuit of length i ending in gate G */
    for (it = cliff_temp[i-1].begin(); it != cliff_temp[i-1].end(); it++) {
      flg = false;
      if (it->second != NULL) {
        flg = nontrivial_id(G, it->second);
      }
      if (!flg) {
        /* Create a circuit for the gate */
        tmp_circ = new Circuit;
        tmp_circ->G = G;
        tmp_circ->next = it->second;
        flg = insert_tree(tmp_circ, cliff_temp, i, false);
        ret = ret || flg;
      }
    }
  }

  return ret;
}

/* Generate the iteration set */
void generate_base_circuits(bool cliffords) {
  Circuit * ins;
  int j;
  Gate G;

  cout << "Generating iteration set V(n, G)\n";
  if (!cliffords) {
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
      cout << "Cliffords of length " << i << "\n";
      flg = generate_cliff(i);
      cout << cliff_temp[i].size() << "\n";
      i++;
    }

    if (i >= CLIFF) cout << "ERROR: unique cliffords of length > " << CLIFF << "\n";
    for (j = 1; j < i; j++) {
      for (k = 0; k < pow(2, num_qubits); k++) {
        for (l = 0; l < num_qubits; l++) {
          if ((k / (int)pow(2, l)) % 2 == 0) {
            G[l] = I;
          } else {
            G[l] = T;
          }
        }
        for (it = cliff_temp[j].begin(); it != cliff_temp[j].end(); it++) {
          ins = new Circuit(G, it->second);
          ins->next = it->second;
          base_list.push_back(ins);
        }
      }
    }
    delete [] cliff_temp;
  }
  cout << "Iteration set size: " << base_list.size() << "\n";
  
}

int max (int a, int b) {
  if (a > b) return a;
  else return b;
}

void generate_sequences(int i, circuit_list &L) {
  bool flg;
  Circuit * tmp_circ;
  map_iter it;
  circuit_iter c;
  struct triple * trip;
  time_t start, end;
  Gate G;

  if (i == 0) {
    Rmatrix V(dim, dim);

    for (int j = 0; j < num_qubits; j++) {
      G[j] = I;
    }

    G.to_Rmatrix(V);
    circuit_table[i].insert(map_elt(Hash_Rmatrix(V), NULL));
    return;
  }

  cout << "--------------------------------------\n";
  cout << "Generating sequences of length " << i << "\n" << flush;
  time(&start);

  /* Generate all the sequences of length i */
  for (c = L.begin(); c != L.end(); c++) {
    G = (*c)->last();
    /* For each circuit of length i ending in gate G */
    for (it = circuit_table[i-1].begin(); it != circuit_table[i-1].end(); it++) {
      flg = false;
      if (it->second != NULL) {
        flg = nontrivial_id(G, it->second);
      }
      if (!flg) {
        tmp_circ = (*c)->append(it->second);
        insert_tree(tmp_circ, circuit_table, i, SYMMS);
        if (it->second != NULL) {
          tmp_circ = (it->second)->append(*c);
          insert_tree(tmp_circ, circuit_table, i, SYMMS);
        }
      }
    }
  }
  time(&end);
  cout << "Time: " << difftime(end, start) << "s\n";
  cout << "# new unitaries: " << circuit_table[i].size() << "\n";
  cout << "--------------------------------------\n" << flush;
}

bool check_it(Circuit * x, Circuit * y, Rmatrix & target) {
  if (dim == reduced_dim) return true;

  Rmatrix left(dim, dim), right(dim, dim);
  x->to_Rmatrix(left);
  y = y->adj(NULL);
  y->to_Rmatrix(right);

  Rmatrix tmp1(dim, dim), tmp2(reduced_dim, reduced_dim);
  tmp1 = left * right;
  tmp1.submatrix(0, 0, reduced_dim, reduced_dim, tmp2);
  return target.phase_eq(tmp2);

}

void exact_search(Rmatrix & U, circuit_list &L) {
  int i, j, k;
  Circuit * ans, * tmp_circ;
  Rmatrix V(dim, dim);
  Rmatrix W(reduced_dim, dim);
  map_iter it;
  Canon canon_form1, canon_form2;
  struct triple * trip;
  int s;
  map_t * mp = (dim == reduced_dim) ? circuit_table : left_table;

  if (SYMMS) s = 2*num_swaps;
  else s = 1;

  time_t start, end;

  generate_sequences(0, L);
  for (i = 1; i < MAX_SEQ; i++) {
    generate_sequences(i, L);
    cout << "Looking for circuits...\n";
    time(&start);
    /* Meet in the middle - Sequences of length 2i + {0, 1} */
    for (j = max(i-1, 1); j <= i; j++) {
      for (it = circuit_table[j].begin(); it != circuit_table[j].end(); it++) {
        for (k = 0; k < s; k++) {

          if (k == 0) {
            tmp_circ = it->second;
          } else if (k % 2 == 0) {
            /* Adjoint case */
            tmp_circ = (it->second)->permute(k/2);
          } else {
            /* Non-adjoint case */
            tmp_circ = (it->second)->permute_adj(k/2, NULL);
          }

          tmp_circ->to_Rmatrix(V, true);
          if (dim != reduced_dim) {
            V.submatrix(0, 0, reduced_dim, dim, W);
            canon_form1 = canonicalize(U*W, SYMMS);
            canon_form2 = canonicalize(W*U, SYMMS);
          } else {
            canon_form1 = canonicalize(U*V, SYMMS);
            canon_form2 = canonicalize(V*U, SYMMS);
          }

          while(!canon_form1.empty()) {
            trip = &(canon_form1.front());
            ans = find_unitary(trip->key, trip->mat, mp[i]).second;
            if (ans != NULL && check_it(ans, tmp_circ, U)) {
              ans->print(tmp_circ);
              cout << "\n" << flush;
            }
            canon_form1.pop_front();
          }
          while(!canon_form2.empty()) {
            trip = &(canon_form2.front());
            ans = find_unitary(trip->key, trip->mat, mp[i]).second;
            if (ans != NULL && check_it(ans, tmp_circ, U)) {
              tmp_circ->print(ans);
              cout << "\n" << flush;
            }
            canon_form2.pop_front();
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

void find_first_nlr() {
  int i, j, k;
  Circuit * ans, * tmp_circ;
  map_iter it, ti;
  Canon canon_form;
  struct triple * trip;
  int s;

  init(3, 3);
  generate_base_circuits(false);
  Rmatrix V(dim, dim);

  // Testing
  Circuit * tst = new Circuit;
  tst->G[0] = X;
  tst->G[1] = X;
  tst->G[2] = X;
  tst->to_Rmatrix(V);
  assert(!V.is_nonlinear_reversible());
  //////////

  for (i = 0; i < MAX_SEQ; i++) {
    generate_sequences(i, base_list);
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

void Rz() {
  init(1, 1);
  generate_base_circuits(false);
  numcorrect = 0;
  numcollision = 0;

  Circuit *h1 = new Circuit(H, new Circuit(T, new Circuit(T, new Circuit(T, NULL))));
  Circuit *h2 = new Circuit(H, new Circuit(Td, new Circuit(H, new Circuit(Td, new Circuit(H, new Circuit(T, h1))))));
  Circuit *h3 = new Circuit(H, new Circuit(T, new Circuit(T, new Circuit(T, h2))));
  Circuit *h4 = new Circuit(H, new Circuit(Td, new Circuit(Td, new Circuit(Td, h3))));
  Circuit *h5 = new Circuit(H, new Circuit(T, new Circuit(T, new Circuit(T, h4))));
  Circuit *h6 = new Circuit(H, new Circuit(T, h5));
  Circuit *h7 = new Circuit(H, new Circuit(T, new Circuit(T, new Circuit(T, h6))));
  Circuit *h8 = new Circuit(H, new Circuit(T, new Circuit(T, new Circuit(T, h7))));
  Circuit *h9 = new Circuit(H, new Circuit(T, new Circuit(T, new Circuit(T, h8))));
  Circuit *fn = new Circuit(X, new Circuit(T, h9));
  Circuit *fin = fn;

  Rmatrix U(dim, dim);
  Unitary V(dim, dim);
  fin->to_Rmatrix(U);
  U.to_Unitary(V);
  exact_search(U, base_list);
}

void Tof() {
  init(3, 3);
  generate_base_circuits(false);
  numcorrect = 0;
  numcollision = 0;

  Circuit * x = new Circuit;
  x->G[0] = C(2);
  x->G[1] = C(2); 
  x->G[2] = X; 

  Rmatrix U(dim, dim);
  x->to_Rmatrix(U);
  exact_search(U, base_list);
}

void CH() {
  init(2, 2);
  generate_base_circuits(false);
  numcorrect = 0;
  numcollision = 0;

  Circuit * x = new Circuit;
  x->G[0] = C(1);
  x->G[1] = H; 

  Rmatrix U(dim, dim);
  x->to_Rmatrix(U);
  exact_search(U, base_list);
}

void TST() {
  init(3, 2);
  generate_base_circuits(false);
  numcorrect = 0;
  numcollision = 0;

  Circuit * a = new Circuit;
  a->G[0] = I;
  a->G[1] = I;
  a->G[2] = X;

  Rmatrix U(dim, dim);
  Rmatrix V(reduced_dim, reduced_dim);
  a->to_Rmatrix(U);
  U.submatrix(0, 0, reduced_dim, reduced_dim, V);
  exact_search(V, base_list);
}

void QFT() {
  init(3, 3);
  generate_base_circuits(false);
  numcorrect = 0;
  numcollision = 0;

  Circuit * a = new Circuit;
  Circuit * b = new Circuit;
  Circuit * c = new Circuit;
  Circuit * d = new Circuit;
  Circuit * e = new Circuit;
  Circuit * f = new Circuit;

  a->G[0] = I;
  a->G[1] = I;
  a->G[2] = H;
  b->G[0] = I;
  b->G[1] = C(2);
  b->G[2] = Td;
  c->G[0] = C(2);
  c->G[1] = I;
  c->G[2] = Sd;
  d->G[0] = I;
  d->G[1] = H;
  d->G[2] = I;
  e->G[0] = C(1);
  e->G[1] = Sd;
  e->G[2] = I;
  f->G[0] = H;
  f->G[1] = I;
  f->G[2] = I;

  f->next = e;
  e->next = d;
  d->next = c;
  c->next = b;
  b->next = a;

  Rmatrix U(dim, dim);
  f->to_Rmatrix(U);
  f->print();
  U.print();
  exact_search(U, base_list);
}

void mem_test() {
  circuit_iter c;
  init(3, 3);
  generate_base_circuits(true);
//  generate_sequences(0, base_list);
 // generate_sequences(1, base_list);
  for (c = base_list.begin(); c != base_list.end(); c++) {
    delete (*c);
  }
}

void test_all() {
  init(3, 3);
  ring_test();
  matrix_test();
  gate_test();
  circuit_test();
}

int main() {
  Tof();

  return 0;
}
