#include "circuit.h"
#include <string>
#include <pthread.h>

#define NUM_THREADS 4

#if ORDERED
typedef multimap<hash_t, Circuit *, cmp_hash> map_t;
typedef map_t::iterator                       map_iter;
typedef map_t::const_iterator           const_map_iter;
#else
#include <unordered_map>
typedef unordered_multimap<hash_t, 
                           Circuit *,
                           hasher,
                           eq_hash>           map_t;
typedef map_t::iterator                       map_iter;
#endif
typedef pair    <hash_t, Circuit *>           map_elt;

typedef list<Circuit *>        circuit_list;
typedef circuit_list::iterator circuit_iter;

typedef struct result {
  bool first;
  Circuit * second;
  map_iter third;
  result() { }
  result(bool b, Circuit * c) { first = b; second = c; }
  result(bool b, Circuit * c, map_iter it) { first = b; second = c; third = it; }
} result;

/* List of base circuits */
circuit_list base_list;

/* Store all unique unitaries up to MAX_SEQ */
map_t circuit_table[MAX_SEQ];
/* Store all unique unitaries with key given by their mapping to the |0> subspace */
map_t left_table[MAX_SEQ];
map_t * cliff_temp;

unsigned int numcorrect = 0, numcollision = 0, numsearch = 0;
unsigned int reserve_num[MAX_SEQ] = {1, 32, 991, 40896, 1418930};

/*-------- threading stuff */
pthread_cond_t data_ready;
pthread_cond_t thrd_ready;
pthread_mutex_t data_lock;
pthread_mutex_t prnt_lock;

Circuit * data_circ = NULL;
int       data_k = 0;
map_t   * data_map = NULL;
bool      data_lnrcosets = false;
bool      data_avail = false;
/*-------------------------------------*/

/*-------------------------------------*/
void output_map_elt(ofstream & out, map_iter t, int depth) {
  output_key(out, t->first);
  (t->second)->output(out, depth);
}

map_elt input_map_elt(ifstream & in, int depth) {
  hash_t key;
  Circuit * circ = new Circuit;
  input_key(in, key);
  circ->input(in, depth);
  return map_elt(key, circ);
}

void output_map(ofstream & out, map_t & map, int depth) {
  map_iter it;
  for (it = map.begin(); it != map.end(); ++it) {
    output_map_elt(out, it, depth);
  }
}

void input_map(ifstream & in, map_t & map, int depth) {
  map_iter it = map.begin();
  while(!(in.peek() == std::ifstream::traits_type::eof())) {
    it = map.insert(it, input_map_elt(in, depth));
  }
}

string gen_filename(int q, int p, int d) {
  stringstream ret;
  ret << "data" << q << "q" << p << "p" << d << "d";
  return ret.str();
}

string gen_cliff_filename(int q, int d) {
  stringstream ret;
  ret << "clifford" << q << "q" << d << "d";
  return ret.str();
}
/*-------------------------------------*/

/* Determine if there is a nontrivial pair of gates that multiply to the identity */
bool nontrivial_id(const Gate & G, const Circuit * circ) {
  int i, j, mask = ~((int)0), tmp = 0, x;
  bool ret = false;
  for (i = 0; i < num_qubits && !ret; i++) {
    if ((1 << i) & mask) {
      if (IS_C(circ->G[i])) {
        x = ~(1 << GET_TARGET(circ->G[i]));
        mask &= x;
        tmp &= x; 
        ret = ret || (G[i] == circ->G[i]);
      } else if (IS_C(G[i])) {
        x = ~(1 << GET_TARGET(G[i]));
        mask &= x;
        tmp &= x; 
        ret = ret || (G[i] == circ->G[i]);
      } else if (G[i] == X && circ->G[i] == X) {
        tmp |= 1 << i;
      } else {
        ret = ret || (G[i] == adjoint[circ->G[i]]);
      }
    }
    /*
    if (!IS_C(circ->G[i])) {
      if (G[i] == X && circ->G[i] == X) {
        ret = true;
        for (j = 0; j < num_qubits; j++) {
          ret = ret && !((G[j] == C(i)) xor (circ->G[j] == C(i)));
        }
      } else { 
        ret = ret || (G[i] == adjoint[circ->G[i]]);
      }
    }
    */
  }
  return ret || (bool)tmp;
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
      return result(true, it->second, ret.first);
    }
  }

  return result(false, NULL, ret.first);
}

result find_unitary(const hash_t & key, const Rmatrix & U, map_t &map) {
  pair<map_iter, map_iter> ret;
  map_iter it;
  Rmatrix V(dim, dim), W(U.rows(), U.cols()), *tmp;

  /* Find the range of elements with equal keys */
  //numsearch++;
  ret = map.equal_range(key);
  for (it = ret.first; it != ret.second; ++it) {
    //numcollision++;
    if (CHECK_EQUIV || U.rows() != dim || U.cols() != dim) {
      (it->second)->to_Rmatrix(V);
      // Determine what submatrix we're comparing 
      if (U.rows() != dim || U.cols() != dim) {
        // In general, if U has different dimension and we found it in map, the keys
        // in map should refer to this submatrix of the actual circuit store
        V.submatrix(0, 0, U.rows(), U.cols(), W);
        tmp = &W;
      } else {
        tmp = &V;
      }
      if (U.phase_eq(*tmp)) {
       // numcorrect++;
        return result(true, it->second, ret.first);
      }
    } else {
     // numcorrect++;
      return result(true, it->second, ret.first);
    }
  }

  return result(false, NULL, ret.first);
}

const result find_unitary(const hash_t & key, const Rmatrix & U, const map_t &map) {
  pair<const_map_iter, const_map_iter> ret;
  const_map_iter it;
  Rmatrix V(dim, dim), W(U.rows(), U.cols()), *tmp;

  /* Find the range of elements with equal keys */
  //numsearch++;
  ret = map.equal_range(key);
  for (it = ret.first; it != ret.second; ++it) {
    //numcollision++;
    if (CHECK_EQUIV || U.rows() != dim || U.cols() != dim) {
      (it->second)->to_Rmatrix(V);
      // Determine what submatrix we're comparing 
      if (U.rows() != dim || U.cols() != dim) {
        // In general, if U has different dimension and we found it in map, the keys
        // in map should refer to this submatrix of the actual circuit store
        V.submatrix(0, 0, U.rows(), U.cols(), W);
        tmp = &W;
      } else {
        tmp = &V;
      }
      if (U == *tmp) {
       // numcorrect++;
        return result(true, it->second);
      }
    } else {
     // numcorrect++;
      return result(true, it->second);
    }
  }

  return result(false, NULL);
}

/* Insert a circuit in a forest of trees, up to a specific depth */
bool insert_tree(Circuit * circ, map_t * tree_list, int depth, bool symm) {
  Rmatrix V(dim, dim);
  Canon * canon_form;
  struct triple * trip;
  bool flg = false, ret = false, del = false;
  int i;
  Circuit * ins_circ;
  result res;

  /* Compute the matrix for circuit C */
  circ->to_Rmatrix(V);
  /* Compute the canonical form */
  canon_form = canonicalize(V, symm);

  // Check to see if it's already synthesized
  // We only need to check the canonical form to get all permutations, inversions and phases
  // Also, we store all canonical forms so that when we search, we only have to look at one canonical form
  while (!canon_form->empty()) {
    trip = &(canon_form->front());
    // Search in each tree up to depth
    flg = false;
    for (i = 1; i <= depth && !flg; i++) {
      res = find_unitary(trip->key, trip->mat, tree_list[i]);
      flg = flg || res.first;
    }
    // If none of the searches turned up anything...
    // Our iterator should refer to the depth i tree
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
      tree_list[depth].insert(res.third, map_elt(trip->key, ins_circ));
      ret = true;

      canon_form->pop_front();
    } else {
      canon_form->clear();
    }
  }
  delete canon_form;

  if (!del && symm) delete_circuit(circ);
  return ret;
}

bool insert_tree(const Rmatrix & V, Circuit * circ, map_t * tree_list, int depth, bool symm) {
  Canon * canon_form;
  struct triple * trip;
  bool flg = false, ret = false, del = false;
  int i;
  Circuit * ins_circ;
  result res;

  /* Compute the canonical form */
  canon_form = canonicalize(V, symm);

  // Check to see if it's already synthesized
  // We only need to check the canonical form to get all permutations, inversions and phases
  // Also, we store all canonical forms so that when we search, we only have to look at one canonical form
  while (!canon_form->empty()) {
    trip = &(canon_form->front());
    // Search in each tree up to depth
    flg = false;
    for (i = 1; i <= depth && !flg; i++) {
      res = find_unitary(trip->key, trip->mat, tree_list[i]);
      flg = flg || res.first;
    }
    // If none of the searches turned up anything...
    // Our iterator should refer to the depth i tree
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
      tree_list[depth].insert(res.third, map_elt(trip->key, ins_circ));
      ret = true;

      canon_form->clear();
    } else {
      canon_form->clear();
    }
  }
  delete canon_form;

  if (!del) delete circ;
  return ret;
}

void generate_proj(const map_t & mp, map_t & proj) {
  Rmatrix V(dim, dim), W(reduced_dim, dim);
  const_map_iter it;
  Canon * canon_form;
  struct triple * trip;
  for (it = mp.begin(); it != mp.end(); it++) {
    (it->second)->to_Rmatrix(V);
    V.submatrix(0, 0, reduced_dim, dim, W);
    canon_form = canonicalize(W, SYMMS);
    trip = &(canon_form->front());
    proj.insert(map_elt(trip->key, it->second));
    canon_form->clear();
    delete canon_form;
  }
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
        G.to_Rmatrix(V);
#if SUBSPACE_ABS
        flg = insert_tree(V * it->first, tmp_circ, cliff_temp, i, false);
#else
        flg = insert_tree(tmp_circ, cliff_temp, i, false);
#endif
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
  if (!TDEPTH && !cliffords) {
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
      if (i == 0) {
        flg = generate_cliff(i);
      } else {
        string s = gen_cliff_filename(num_qubits, i);
        ifstream in;
        in.open(s.c_str(), ifstream::binary);
        if (in.peek() != std::ifstream::traits_type::eof()) {
          input_map(in, cliff_temp[i], i);
        } else {
          flg = generate_cliff(i);
          if (flg) {
            ofstream out;
            out.open(s.c_str(), ofstream::binary);
            output_map(out, cliff_temp[i], i);
            out.close();
          }
        }
        in.close();
      }
      cout << cliff_temp[i++].size() << "\n" << flush;
    }

    if (i >= CLIFF) cout << "ERROR: unique cliffords of length > " << CLIFF << "\n";
    for (j = 1; j < i; j++) {
      for (k = 0; k < (1 << num_qubits); k++) {
        for (l = 0; l < num_qubits; l++) {
          if ((k / (1 << l)) % 2 == 0) {
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

void generate_sequences(int i, circuit_list &L) {
  bool flg;
  Circuit * tmp_circ;
  map_iter it;
  circuit_iter c;
  struct triple * trip;
  time_t start, end;
  Gate G;
  Rmatrix V(dim, dim);
  hash_t key;
  string s;

  if (i == 0) {
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

#if !ORDERED
  circuit_table[i].reserve(reserve_num[i]);
#endif
  if (SERIALIZE && i != 0) {
    s = gen_filename(num_qubits, num_qubits, i);
    ifstream in;
    in.open(s.c_str(), ifstream::binary);
    if (in.peek() != std::ifstream::traits_type::eof()) {
      input_map(in, circuit_table[i], i);
    } else {
      ofstream out;
      out.open(s.c_str(), ofstream::binary);

      /* Generate all the sequences of length i */
      for (c = L.begin(); c != L.end(); c++) {
        /* For each circuit of length i ending in gate G */
        for (it = circuit_table[i-1].begin(); it != circuit_table[i-1].end(); it++) {
          flg = false;
          if (it->second != NULL) {
            flg = nontrivial_id((*c)->last(), it->second);
          }
          if (!flg) {
            tmp_circ = (*c)->append(it->second);
#if SUBSPACE_ABS
            (*c)->to_Rmatrix(V);
            key = Hash_Rmatrix(V);
            insert_tree(key * it->first, tmp_circ, circuit_table, i, SYMMS);
#else
            insert_tree(tmp_circ, circuit_table, i, SYMMS);
#endif
            if (it->second != NULL) {
              tmp_circ = (it->second)->append(*c);
#if SUBSPACE_ABS
              insert_tree(it->first * key, tmp_circ, circuit_table, i, SYMMS);
#else
              insert_tree(tmp_circ, circuit_table, i, SYMMS);
#endif
            }
          }
        }
      }

      output_map(out, circuit_table[i], i);
      out.close();
    }
    in.close();
  } else {
    /* Generate all the sequences of length i */
    for (c = L.begin(); c != L.end(); c++) {
      /* For each circuit of length i ending in gate G */
      for (it = circuit_table[i-1].begin(); it != circuit_table[i-1].end(); it++) {
        flg = false;
        if (it->second != NULL) {
          flg = nontrivial_id((*c)->last(), it->second);
        }
        if (!flg) {
          tmp_circ = (*c)->append(it->second);
#if SUBSPACE_ABS
          (*c)->to_Rmatrix(V);
          key = Hash_Rmatrix(V);
          insert_tree(key * it->first, tmp_circ, circuit_table, i, SYMMS);
#else
          insert_tree(tmp_circ, circuit_table, i, SYMMS);
#endif
          if (it->second != NULL) {
            tmp_circ = (it->second)->append(*c);
#if SUBSPACE_ABS
            insert_tree(it->first * key, tmp_circ, circuit_table, i, SYMMS);
#else
            insert_tree(tmp_circ, circuit_table, i, SYMMS);
#endif
          }
        }
      }
    }
  }

  if (i != 0 && num_qubits != num_qubits_proj) {
    if (SERIALIZE) {
      s = gen_filename(num_qubits, num_qubits_proj, i);
      ifstream in;
      in.open(s.c_str(), ifstream::binary);
      if (in.peek() != std::ifstream::traits_type::eof()) {
        input_map(in, left_table[i], i);
      } else {
        generate_proj(circuit_table[i], left_table[i]);
        ofstream out;
        out.open(s.c_str(), ofstream::binary);
        output_map(out, left_table[i], i);
        out.close();
      }
      in.close();
    } else {
      generate_proj(circuit_table[i], left_table[i]);
    }
  }

  time(&end);
  cout << "Time: " << difftime(end, start) << "s\n";
  cout << "# new unitaries: " << circuit_table[i].size() << "\n";
  cout << "# searches so far: " << numsearch << "\n";
  cout << "equivalent unitary vs equivalent key: " << numcorrect << " / " << numcollision << "\n";
#if !ORDERED
  int m = 0;
  for (int j = 0; j < circuit_table[i].bucket_count(); j++) {
    m = max(m, circuit_table[i].bucket_size(j));
  }
  cout << "Bucket count: " << circuit_table[i].bucket_count() << "\n";
  cout << "Load factor: " << circuit_table[i].load_factor() << "\n";
  cout << "Max bucket size: " << m << "\n";
#endif
  cout << "--------------------------------------\n" << flush;
}

bool check_it(const Circuit * x, const Circuit * y, const Rmatrix & target) {
  Rmatrix tmp1(dim, dim), tmp2(dim, dim);
  x->to_Rmatrix(tmp1);
  y->to_Rmatrix(tmp2);

  tmp1 *= tmp2;
  if (dim == reduced_dim) {
    if (!target.phase_eq(tmp1)) {
      cout << "WARNING: may not be correct\n" << flush;
    }
    return true;
  } else {
    tmp1.submatrix(0, 0, reduced_dim, reduced_dim, tmp2);
    return target.phase_eq(tmp2);
  }

}

void worker_thrd(const Rmatrix * arg) {
  Circuit * circ, * tmp_circ, *tmp_circ2;
  const Circuit * ans;
  int k, i, j;
  bool lnrcosets;
  map_t * mp;
  Rmatrix V(dim, dim), W(dim, dim), U = *arg;
  Canon * canon_form1, * canon_form2;
  struct triple * trip;

  pthread_mutex_lock(&data_lock);
  while(1) {
    if (data_avail) {
      // Copy our data
      circ = data_circ;
      k = data_k;
      mp = data_map;
      lnrcosets = data_lnrcosets;
      data_avail = false;
      // Signal the mamma thread
      pthread_cond_signal(&thrd_ready);
      // Unlock the lock
      pthread_mutex_unlock(&data_lock);

      // Perform the search
      if (k == 0) {
        tmp_circ = circ;
      } else if (k == 1) {
        tmp_circ = (circ)->adj(NULL);
      } else if (k % 2 == 0) {
        tmp_circ = (circ)->permute(k/2);
      } else {
        tmp_circ = (circ)->permute_adj(k/2, NULL);
      }

      tmp_circ->to_Rmatrix(V, true);
      if (dim != reduced_dim) {
        V.submatrix(0, 0, reduced_dim, dim, W);
        canon_form1 = canonicalize(U*W, true);
      } else {
        canon_form1 = canonicalize(U*V, SYMMS);
      }

      trip = &(canon_form1->front());
      ans = find_unitary(trip->key, trip->mat, *mp).second;
      if (ans != NULL) {
        if (trip->adjoint == false) {
          tmp_circ2 = ans->permute(trip->permutation, true);
        } else {
          tmp_circ2 = ans->permute_adj(trip->permutation, true, NULL);
        }
        if (check_it(tmp_circ2, tmp_circ, U)) {
          pthread_mutex_lock(&prnt_lock);
          tmp_circ2->print(tmp_circ);
          cout << "\n" << flush;
          pthread_mutex_unlock(&prnt_lock);
        }
        delete_circuit(tmp_circ2);
      }
      canon_form1->clear();
      delete canon_form1;

      if (lnrcosets && num_qubits == num_qubits_proj) {
        if (dim != reduced_dim) {
          canon_form2 = canonicalize(W*U, SYMMS);
        } else {
          canon_form2 = canonicalize(V*U, SYMMS);
        }

        trip = &(canon_form2->front());
        ans = find_unitary(trip->key, trip->mat, *mp).second;
        if (ans != NULL) {
          if (trip->adjoint == false) {
            tmp_circ2 = ans->permute(trip->permutation, true);
          } else {
            tmp_circ2 = ans->permute_adj(trip->permutation, true, NULL);
          }
          if (check_it(tmp_circ, tmp_circ2, U)) {
            pthread_mutex_lock(&prnt_lock);
            tmp_circ->print(tmp_circ2);
            cout << "\n" << flush;
            pthread_mutex_unlock(&prnt_lock);
          }
          delete_circuit(tmp_circ2);
        }
        canon_form2->clear();
        delete canon_form2;
      }
      if (k != 0) delete_circuit(tmp_circ);

      // Lock before we reenter the loop
      pthread_mutex_lock(&data_lock);
    } else {
      // Wait until data is ready
      pthread_cond_wait(&data_ready, &data_lock);
    }
  }
  pthread_exit(NULL);
}

void exact_search(Rmatrix & U, circuit_list &L) {
  int i, j, k;
  Circuit * ans, * tmp_circ, * tmp_circ2;
  Rmatrix V(dim, dim);
  Rmatrix W(reduced_dim, dim);
  map_iter it;
  Canon canon_form1, canon_form2;
  struct triple * trip;
  int s;
  map_t * mp = (dim == reduced_dim) ? circuit_table : left_table;

  if (SYMMS) s = 2*num_swaps;
  else s = 1;

#if SUBSPACE_ABS
  hash_t Uhash = Hash_Rmatrix(U);
#endif

  time_t start, end;

  //Threading stuff
  pthread_mutex_init(&data_lock, NULL);
  pthread_mutex_init(&prnt_lock, NULL);
  pthread_cond_init(&data_ready, NULL);
  pthread_cond_init(&thrd_ready, NULL);

  pthread_t thrds[NUM_THREADS];
  for (i = 0; i < NUM_THREADS; i++) {
    pthread_create(thrds + i, NULL, &(worker_thrd), (void *)(&U));
  }
  //---------------------

  generate_sequences(0, L);
  pthread_mutex_lock(&data_lock);
  for (i = 1; i < MAX_SEQ; i++) {
    generate_sequences(i, L);
    cout << "Looking for circuits...\n";
    time(&start);
    /* Meet in the middle - Sequences of length 2i + {0, 1} */
    for (j = max(i-1, 1); j <= i; j++) {
      for (it = circuit_table[j].begin(); it != circuit_table[j].end(); it++) {
        for (k = 0; k < s; k++) {
          // Wait until we can write data
          if (!data_avail) {
            // Write data
            data_circ = it->second;
            data_k = k;
            data_lnrcosets = (i != j);
            data_map = mp + i;
            data_avail = true;
            // Signal workers that data is ready
            pthread_cond_signal(&data_ready);
          }
          pthread_cond_wait(&thrd_ready, &data_lock);
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
  init_ht();
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
  init_ht();
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
  init_ht();
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
  generate_base_circuits(false);
  generate_sequences(0, base_list);
  generate_sequences(1, base_list);
  generate_sequences(2, base_list);
  generate_sequences(3, base_list);
  /*
  for (c = base_list.begin(); c != base_list.end(); c++) {
    delete (*c);
  }
  */
  int i = 1;
  while(1) {
    i++;
  }
}

void test_all() {
  init(3, 3);
  ring_test();
  matrix_test();
  gate_test();
  circuit_test();
}

void mem_tst() {
  init(3, 3);
  init_ht();
  generate_base_circuits(false);
}

int main() {
  Tof();

  return 0;
}
