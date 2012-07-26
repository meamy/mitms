#include "search.h"
#include <string>
#include <pthread.h>
#include <iomanip>

/* Circuit databases */
typedef multimap<hash_t, Circuit, cmp_hash> map_t;
typedef map_t::iterator                     map_iter;
typedef map_t::const_iterator         const_map_iter;
typedef pair    <hash_t, Circuit>           map_elt;

/* Circuit lists */
typedef list<Circuit>          circuit_list;
typedef circuit_list::iterator circuit_iter;

/* Ordered Circuit lists */
typedef multimap<int, Circuit>     ord_circuit_list;
typedef ord_circuit_list::iterator ord_circuit_iter;
typedef pair    <int, Circuit>     ord_circuit_pair;

/* Structure for returning results of searches */
typedef struct result {
  bool first;
  Circuit second;
  map_iter third;
  result() { }
  result(bool b, Circuit c) { first = b; second = c; }
  result(bool b, Circuit c, map_iter it) { first = b; second = c; third = it; }
} result;

/* List of the clifford group circuits */
circuit_list * cliff_list;

unsigned int numcorrect = 0, numcollision = 0, numsearch = 0;

/*-------- threading stuff */
pthread_cond_t data_ready;
pthread_cond_t thrd_ready;
pthread_mutex_t data_lock;
pthread_mutex_t prnt_lock;
pthread_mutex_t map_lock;

map_t            * data_map = NULL;
ord_circuit_list * data_res = NULL;
Circuit            data_circ;
Rmatrix            data_mat;
int                data_k = 0;
bool               data_left = false;
bool               data_avail = false;
int                data_num = 0;
bool               data_exit = false;
/*-------------------------------------*/


/* It seems like this function could be expanded to deal with things like phase, 
   but DO NOT do so -- this way is more efficient */
result find_unitary(const hash_t & key, const Unitary & U, map_t &map) {
  pair<map_iter, map_iter> ret;
  map_iter it, ti;
  Unitary V(dim, dim);

  ret = map.equal_range(key);
  for (it = ret.first; it != ret.second; ++it) {
    (it->second).to_Unitary(V);
    if (dist(U, V) < config::precision) {
      return result(true, it->second, ret.first);
    }
  }

  return result(false, Circuit(), ret.first);
}

result find_unitary(const hash_t & key, const Rmatrix & U, map_t & map) {
  pair<map_iter, map_iter> ret;
  map_iter it;
  Rmatrix V(dim, dim), W(U.rows(), U.cols()), *tmp;

  /* Find the range of elements with equal keys */
  numsearch++;
  ret = map.equal_range(key);
  for (it = ret.first; it != ret.second; it++) {
    numcollision++;
    if (config::check_equiv || U.rows() != dim || U.cols() != dim) {
      (it->second).to_Rmatrix(V);
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
        numcorrect++;
        return result(true, it->second, it);
      }
    } else {
      numcorrect++;
      return result(true, it->second, it);
    }
  }

  return result(false, Circuit(), ret.first);
}

inline const result find_unitary(const hash_t & key, const Rmatrix & U, const map_t & map) {
  return find_unitary(key, U, map);
}

/*-------------------------------------*/
void output_map_elt(ofstream & out, map_iter t, int depth) {
  output_key(out, t->first);
  (t->second).output(out, depth);
}

map_elt input_map_elt(ifstream & in, int depth) {
  hash_t key;
  Circuit circ;
  input_key(in, key);
  circ.input(in, depth);
  return map_elt(key, circ);
}

void output_map(ofstream & out, map_t * map, int depth) {
  map_iter it;
  map_iter lst = map[depth].begin();
  Rmatrix U(dim, dim), V(dim, dim);

  for (it = map[depth].begin(); it != map[depth].end(); it++) {
    if (lst == it || (!(lst->first == it->first))) {
      output_map_elt(out, it, depth);
      lst = it;
    } else {
      lst->second.to_Rmatrix(U);
      it->second.to_Rmatrix(V);
      if (!U.phase_eq(V)) {
        output_map_elt(out, it, depth);
        lst = it;
      } else {
        lst->second.print();
        U.print();
        it->second.print();
        V.print();
        cout << "\n";
      }
    }
  }
}

void input_map(ifstream & in, map_t * map, int depth) {
  map_iter it = map[depth].begin();
  hash_t tmp_key;
  Rmatrix tmp_mat(dim, dim);
  result res;

  while(!(in.peek() == std::ifstream::traits_type::eof())) {
    it = map[depth].insert(it, input_map_elt(in, depth));
  }
}

/*-------------------------------------*/

/* Insert a circuit in a forest of trees, up to a specific depth */
void insert_tree(const Rmatrix & V, Circuit & circ, map_t * tree_list, int depth, bool perms, bool inv) {
  Canon * canon_form;
  struct triple * trip;
  bool flg = false, del = false;
  int i;
  Circuit ins_circ;
  result res;

  /* Compute the canonical form */
  canon_form = canonicalize(V, config::mod_phase, perms, inv);

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
    if (!flg || ((i > depth) && (circ.cost() < res.second.cost()))) {
      if (flg) {
        tree_list[i-1].erase(res.third);
      }
      if (trip->permutation == 0 && trip->adjoint == false) {
        ins_circ = circ; 
        del = true;
      } else {
        ins_circ = circ.transform(trip->permutation, trip->adjoint);
      }

      if (flg) {
        tree_list[depth].insert(map_elt(trip->key, ins_circ));
      } else {
        tree_list[depth].insert(res.third, map_elt(trip->key, ins_circ));
      }

      canon_form->pop_front();
    } else {

      canon_form->clear();
    }
  }
  delete canon_form;

  if (!del) delete_circuit(circ);
}

inline void insert_tree(Circuit & circ, map_t * tree_list, int depth, bool perms, bool inv) {
  Rmatrix V(dim, dim);
  circ.to_Rmatrix(V);
  insert_tree(V, circ, tree_list, depth, perms, inv);
}

/* Threaded version -- for generating cliffords */
void insert_tree_thrd(void * arguments) {
  Rmatrix V(dim, dim);
  Canon * canon_form;
  struct triple * trip;
  bool flg = false, del = false;
  int i, depth;
  Circuit circ, ins_circ;
  result res;
  map_t * tree_list;
  bool * symms = arguments;

  pthread_mutex_lock(&data_lock);
  while(1) {
    if (data_exit) {
      pthread_mutex_unlock(&data_lock);
      pthread_exit(NULL);
    } else if (data_avail) {
      // Copy our data
      circ = data_circ;
      depth = data_k;
      tree_list = data_map;
      data_avail = false;
      data_num++;
      // Signal the mamma thread
      pthread_cond_signal(&thrd_ready);
      // Unlock the lock
      pthread_mutex_unlock(&data_lock);

      /* Compute the matrix for circuit C */
      circ.to_Rmatrix(V);
      /* Compute the canonical form */
      canon_form = canonicalize(V, config::mod_phase, symms[0], symms[1]);

      // Check to see if it's already synthesized
      // We only need to check the canonical form to get all permutations, inversions and phases
      // Also, we store all canonical forms so that when we search, we only have to look at one canonical form
      del = false;
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
          if (trip->permutation == 0 && trip->adjoint == false) {
            ins_circ = circ; 
            del = true;
          } else {
            ins_circ = circ.transform(trip->permutation, trip->adjoint);
          }

          pthread_mutex_lock(&map_lock);
          tree_list[depth].insert(res.third, map_elt(trip->key, ins_circ));
          pthread_mutex_unlock(&map_lock);

          canon_form->pop_front();
        } else {
          canon_form->clear();
        }
      }
      delete canon_form;

      if (!del) delete_circuit(circ);

      // Lock before we reenter the loop
      pthread_mutex_lock(&data_lock);
      data_num--;
    } else {
      // Wait until data is ready
      pthread_cond_wait(&data_ready, &data_lock);
    }
  }
}

void generate_proj(int depth, map_t * mp, map_t * left_table) {
  Rmatrix V(dim, dim), W(dim_proj, dim);
  map_iter it;
  Circuit tmp_circ;
  int k, num = 0, p = 0;

	cout << "|";

  for (it = mp->begin(); it != mp->end(); it++) {
    for (k = 0; k < 2*num_perms; k++) {
      tmp_circ = (it->second).transform(k/2, k % 2);
      tmp_circ.to_Rmatrix(V);
      V.submatrix(0, 0, dim_proj, dim, W);
      insert_tree(W, tmp_circ, left_table, depth, false, false);
    }
		delete_circuit(it->second);
		num++;
		if (num*37 / mp->size() >= p) {
			cout << "=" << flush;
			p++;
		}
  }
	cout << "|\n";
}

void generate_cliff(int i, map_t * cliff_temp) {
  int j, k;
  Circuit gt(1);
  Rmatrix V(dim, dim);
  map_iter it;
  bool flg = false;

  // Reset gate G
  if (i == 0) {
    gt.to_Rmatrix(V);
    cliff_temp[i].insert(map_elt(Hash_Rmatrix(V), Circuit()));
    return;
  }
  else if (i == 1) {
    gt.to_Rmatrix(V);
    cliff_temp[i].insert(map_elt(Hash_Rmatrix(V), gt));
  }

  /* Generate all the sequences of length i */
  pthread_mutex_lock(&data_lock);
  data_map = cliff_temp;
  while(!((gt.cliffpp()).eye())) {
    /* For each circuit of length i ending in gate G */
    for (it = cliff_temp[i-1].begin(); it != cliff_temp[i-1].end(); it++) {
      // Wait until we can write data
      if (!data_avail) {
        flg = false;
        if (!(it->second.empty())) {
          flg = nontrivial_id(gt.last(), (it->second).first());
        }
        if (!flg) {
          // Write data
          data_circ = gt.append(it->second);
          data_k = i;
          data_avail = true;
          // Signal workers that data is ready
          pthread_cond_signal(&data_ready);
          pthread_cond_wait(&thrd_ready, &data_lock);
        }
      }
    }
  }
  // Wait for all threads to finish 
  while (data_num > 0) {
    pthread_mutex_unlock(&data_lock);
    pthread_mutex_lock(&data_lock);
  }
}

/* Generate the iteration set */
circuit_list * generate_base_circuits() {
  circuit_list * ret = new circuit_list;
  Circuit ins;
  int i, j;
  Circuit gt(1);

  cout << "Generating iteration set V(n, G)\n";
  if (!config::tdepth) {
    // Reset gate G
    ret->push_back(gt.copy());

    /* Generate all the sequences of length i */
    while(!((++gt).eye())) {
      ret->push_back(gt.copy());
    }
  } else {
    //Threading stuff
    pthread_mutex_init(&data_lock, NULL);
    pthread_mutex_init(&prnt_lock, NULL);
    pthread_mutex_init(&map_lock, NULL);
    pthread_cond_init(&data_ready, NULL);
    pthread_cond_init(&thrd_ready, NULL);

    pthread_t * thrds = new pthread_t[config::num_threads];
    bool symms[2] = {false, false};
    for (i = 0; i < config::num_threads; i++) {
      pthread_create(thrds + i, NULL, &(insert_tree_thrd), symms);
    }
    //---------------------
  
    int k, l;
    map_iter it;

    map_t * cliff_temp = new map_t[config::max_seq];
    i = 0;
    generate_cliff(i++, cliff_temp);
    while(cliff_temp[i-1].size() != 0 && i < config::max_seq) {
      cout << "Cliffords of length " << i << "\n";

      string s = gen_cliff_filename(num_qubits, i);
      ifstream in;
      in.open(s.c_str(), ifstream::binary);
      if (in.peek() != std::ifstream::traits_type::eof()) {
        input_map(in, cliff_temp, i);
      } else {
        generate_cliff(i, cliff_temp);
        while (data_num > 0) {
          pthread_mutex_unlock(&data_lock);
          pthread_mutex_lock(&data_lock);
        }
        pthread_mutex_unlock(&data_lock);
        if (cliff_temp[i].size() != 0) {
          ofstream out;
          out.open(s.c_str(), ofstream::binary);
          output_map(out, cliff_temp, i);
          out.close();
        }
      }
      in.close();
      cout << cliff_temp[i++].size() << "\n" << flush;
    }

    pthread_mutex_lock(&data_lock);
    data_exit = true;
    pthread_cond_broadcast(&data_ready);
    pthread_mutex_unlock(&data_lock);
    for (j = 0; j < config::num_threads; j++) {
      pthread_join(thrds[j], NULL);
    }
    if (i >= config::max_seq) cout << "ERROR: unique cliffords of length > " << config::max_seq << "\n";

    /* Generate C_n */
    cliff_list = new circuit_list;
    for (j = 1; j < i; j++) {
      for (it = cliff_temp[j].begin(); it != cliff_temp[j].end(); it++) {
        cliff_list->push_back(it->second);
      }
    }
    cout << "Clifford group size: " << cliff_list->size() << "\n";
    delete [] cliff_temp;

    /* Generate V_n,G */
    gate G = gt.last();
    for (k = 1; k < (1 << num_qubits); k++) {
      for (l = 0; l < num_qubits; l++) {
        if ((k / (1 << l)) % 2 == 0) {
          G[l] = I;
        } else {
          G[l] = T;
        }
      }
      ins = gt.copy();
      ret->push_back(ins);
    }
  }
  cout << "Iteration set size: " << ret->size() << "\n";
  return ret; 
}

void generate_sequences(int i, circuit_list * L, map_t * circ_table) {
  bool flg;
  Circuit tmp_circ;
  map_iter it;
  circuit_iter c;
  Rmatrix V(dim, dim), W(dim, dim);
  int num = 0, p = 0;

	cout << "|";

  /* For each circuit of length i ending in gate G */
  for (it = circ_table[i-1].begin(); it != circ_table[i-1].end(); it++) {
    if (!(it->second.empty())) (it->second).to_Rmatrix(V);
    else V = eye(dim, dim);
    /* Generate all the sequences of length i */
    for (c = L->begin(); c != L->end(); c++) {
      (*c).to_Rmatrix(W);
      flg = false;
      /*
         if (it->second != NULL) {
         flg = nontrivial_id((*c)->last(), (it->second)->G);
         }
       */
      if (!flg) {
        tmp_circ = (*c).append(it->second);
        insert_tree(W*V, tmp_circ, circ_table, i, config::mod_perms, config::mod_invs);
        if (config::mod_invs && !(it->second.empty())) {
          tmp_circ = (it->second).append(*c);
          insert_tree(V*W, tmp_circ, circ_table, i, config::mod_perms, config::mod_invs);
        }
      }
    }
		num++;
		if (num*37 / circ_table[i-1].size() >= p) {
			cout << "=" << flush;
			p++;
		}
  }
	cout << "|\n";
}

void load_sequences(int i, circuit_list * L, map_t * circ_table) {
  struct timespec start, end;
  string s;

  if (i == 0) {
    Circuit gt(1);
    Rmatrix V(dim, dim);

    gt.to_Rmatrix(V);
    circ_table[i].insert(map_elt(Hash_Rmatrix(V), Circuit()));
    return;
  }

  cout << "--------------------------------------\n";
  cout << "Generating sequences of length " << i << "\n" << flush;
  clock_gettime(CLOCK_MONOTONIC, &start);

  if (config::serialize) {
    s = gen_filename(num_qubits, num_qubits, i);
    ifstream in;
    in.open(s.c_str(), ifstream::binary);
    if (in.peek() != std::ifstream::traits_type::eof()) {
      input_map(in, circ_table, i);
    } else {
      ofstream out;
      out.open(s.c_str(), ofstream::binary);

      generate_sequences(i, L, circ_table);

      output_map(out, circ_table, i);
      out.close();
    }
    in.close();
  } else {
    generate_sequences(i, L, circ_table);
  }

  clock_gettime(CLOCK_MONOTONIC, &end);
  cout << fixed << setprecision(3);
  cout << "Time: " << (end.tv_sec + (double)end.tv_nsec/1000000000) - (start.tv_sec + (double)start.tv_nsec/1000000000) << " s\n";
  cout << "# new unitaries: " << circ_table[i].size() << "\n";
  cout << "# searches so far: " << numsearch << "\n";
  cout << "equivalent unitary vs equivalent key: " << numcorrect << " / " << numcollision << "\n";
  cout << "--------------------------------------\n" << flush;
}

void load_proj(int i, circuit_list * L, map_t * circ_table, map_t * left_table) {
  struct timespec start, end;
  string s;

  cout << "--------------------------------------\n";
  cout << "Generating projective sequences of length " << i << "\n" << flush;
  clock_gettime(CLOCK_MONOTONIC, &start);

    if (config::serialize) {
      s = gen_filename(num_qubits, num_qubits_proj, i);
      ifstream in;
      in.open(s.c_str(), ifstream::binary);
      if (in.peek() != std::ifstream::traits_type::eof()) {
        input_map(in, left_table, i);
      } else {
        ofstream out;
        out.open(s.c_str(), ofstream::binary);

				load_sequences(i, L, circ_table);
        generate_proj(i, circ_table + i, left_table);

        output_map(out, left_table, i);
        out.close();
      }
      in.close();
    } else {
			generate_sequences(i, L, circ_table);
      generate_proj(i, circ_table + i, left_table);
    }

  clock_gettime(CLOCK_MONOTONIC, &end);
  cout << fixed << setprecision(3);
  cout << "Time: " << (end.tv_sec + (double)end.tv_nsec/1000000000) - (start.tv_sec + (double)start.tv_nsec/1000000000) << " s\n";
  cout << "# new unitaries: " << left_table[i].size() << "\n";
  cout << "# searches so far: " << numsearch << "\n";
  cout << "equivalent unitary vs equivalent key: " << numcorrect << " / " << numcollision << "\n";
  cout << "--------------------------------------\n" << flush;
}

bool check_it(const Circuit & x, const Circuit & y, const Rmatrix & target) {
  Rmatrix tmp1(dim, dim), tmp2(dim, dim);
  x.to_Rmatrix(tmp1);
  y.to_Rmatrix(tmp2);

  tmp1 *= tmp2;
  if (dim == dim_proj) {
    if (!target.phase_eq(tmp1)) {
      cout << "WARNING: may not be correct\n" << flush;
    }
    return true;
  } else {
    tmp1.submatrix(0, 0, dim_proj, dim_proj, tmp2);
    return target.phase_eq(tmp2);
  }

}

void worker_thrd(const Rmatrix * arg) {
  Circuit circ, tmp_circ, tmp_circ2, tmp_circ3, ans;
  int k, i, j, cst;
  bool left_multiply;
  map_t * mp;
  Rmatrix V(dim, dim), W(dim_proj, dim), U = *arg;
  Canon * canon_form;
  struct triple * trip;

  pthread_mutex_lock(&data_lock);
  while(1) {
    if (data_avail) {
      // Copy our data
      circ = data_circ;
      k = data_k;
      left_multiply = data_left;
      if (k % 2 == 0) {
        data_mat.permute_adj(V, k/2);
      } else {
        data_mat.permute(V, k/2);
      }
      mp = data_map;
      data_avail = false;
      data_num++;
      // Signal the mamma thread
      pthread_cond_signal(&thrd_ready);
      // Unlock the lock
      pthread_mutex_unlock(&data_lock);

      // Perform the search
      if (dim != dim_proj) {
        V.submatrix(0, 0, dim_proj, dim, W);
        canon_form = left_multiply ? canonicalize(W*U, false, false, false) : 
                                     canonicalize(U*W, false, false, false);
      } else {
        canon_form = left_multiply ? canonicalize(V*U) : canonicalize(U*V);
      }

      trip = &(canon_form->front());
      ans = find_unitary(trip->key, trip->mat, *mp).second;
      if (!ans.empty()) {
        // Generate the two circuit halves
        tmp_circ = (k == 0) ? circ : circ.transform(k/2, k % 2 == 1);
        tmp_circ2 = ans.transform(-(trip->permutation), trip->adjoint);

        // Check that the circuit is correct
        if (left_multiply ? check_it(tmp_circ, tmp_circ2, U) : 
                            check_it(tmp_circ2, tmp_circ, U)) {
          tmp_circ3 = left_multiply ? tmp_circ.append(tmp_circ2) : 
                                      tmp_circ2.append(tmp_circ);
          cst = tmp_circ3.cost();
          pthread_mutex_lock(&prnt_lock);
					cout << flush;
          data_res->insert(ord_circuit_pair(cst, tmp_circ3));
          pthread_mutex_unlock(&prnt_lock);
        }

        if (k != 0) delete_circuit(tmp_circ);
        delete_circuit(tmp_circ2);
      }
      canon_form->clear();
      delete canon_form;

      // Lock before we reenter the loop
      pthread_mutex_lock(&data_lock);
      data_num--;
    } else {
      // Wait until data is ready
      pthread_cond_wait(&data_ready, &data_lock);
    }
  }
  pthread_exit(NULL);
}

void exact_search(Rmatrix & U) {
  int num = 0, p = 0;
  int i, j, k;
  int pe = config::mod_perms ? num_perms : 1;
  int in = config::mod_invs  ?         1 : 2;
  map_t * circ_table = new map_t[config::max_seq];
  map_t * left_table = (config::ancilla == 0) ? NULL : new map_t[config::max_seq];
  map_t * mp         = (config::ancilla == 0) ? circ_table : left_table;
  map_iter it;
  circuit_list * base_list;
  ord_circuit_list * res_list = data_res = new ord_circuit_list;
  ord_circuit_iter ti;
  struct timespec start, end;
	data_mat = Rmatrix(dim, dim);

  // Do this first so that the threads don't conflict
  base_list = generate_base_circuits();

  //Initialize the circuit tables
  //Threading stuff
  pthread_mutex_init(&data_lock, NULL);
  pthread_mutex_init(&prnt_lock, NULL);
  pthread_cond_init(&data_ready, NULL);
  pthread_cond_init(&thrd_ready, NULL);

  pthread_t * thrds = new pthread_t[config::num_threads];
  for (i = 0; i < config::num_threads; i++) {
    pthread_create(thrds + i, NULL, &(worker_thrd), (void *)(&U));
  }
  //---------------------

  load_sequences(0, base_list, circ_table);
  pthread_mutex_lock(&data_lock);
  for (i = 1; i < config::max_seq; i++) {
		if (config::ancilla == 0) {
    	load_sequences(i, base_list, circ_table);
		} else {
			load_proj(i, base_list, circ_table, left_table);
		}

    // Meet in the middle - Sequences of length 2i + {0, 1}
    data_map = mp + i;
    for (j = max(i-1, 1); j <= i; j++) {
    	cout << "Looking for circuits with depth " << 2*i - (i - j) << "...\n";
			cout << "|";
			num = 0;
			p = 0;
    	clock_gettime(CLOCK_MONOTONIC, &start);

			// Look for circuits
      for (it = mp[j].begin(); it != mp[j].end(); it++) {
        data_circ = it->second;
        (it->second).to_Rmatrix(data_mat);
        for (k = 0; k < 2*pe; k += in) {
          // Write data
          data_k = k;
          data_avail = true;
          // Signal workers that data is ready
          pthread_cond_signal(&data_ready);
          pthread_cond_wait(&thrd_ready, &data_lock);
        }
				num++;
				if (num*37 / mp[j].size() >= p) {
					cout << "=" << flush;
					p++;
				}
      }
			cout << "|\n";

			// Wait for all threads to finish 
			while (data_num > 0) {
				pthread_mutex_unlock(&data_lock);
				pthread_mutex_lock(&data_lock);
			}
			for (ti = res_list->begin(); ti != res_list->end(); ++ti) {
				(ti->second).print();
				cout << "Cost " << ti->first << "\n\n" << flush;
				delete_circuit(ti->second);
				res_list->erase(ti);
			}

			clock_gettime(CLOCK_MONOTONIC, &end);
			cout << fixed << setprecision(3);
			cout << "Time: " << (end.tv_sec + (double)end.tv_nsec/1000000000) - (start.tv_sec + (double)start.tv_nsec/1000000000) << " s\n";
			cout << "equivalent unitary vs equivalent key: " << numcorrect << " / " << numcollision << "\n";
			cout << "--------------------------------------\n" << flush;
    }
  }

  delete [] circ_table;
  if (left_table != NULL) delete [] left_table;
  delete base_list;
  delete [] thrds;
}

void exact_search_tdepth(Rmatrix & U) {
  int i, j, k;
  int pe = config::mod_perms ? num_perms : 1;
  int in = config::mod_invs  ?         1 : 2;
  map_t * circ_table = new map_t[config::max_seq];
  map_t * left_table = (config::ancilla == 0) ? NULL : new map_t[config::max_seq];
  map_t * mp         = (config::ancilla == 0) ? circ_table : left_table;
  map_iter it;
  circuit_list * base_list;
  ord_circuit_list * res_list = data_res = new ord_circuit_list;
  ord_circuit_iter ti;
  struct timespec start, end;

  circuit_iter c;
  Rmatrix tmp(dim, dim);

  // Do this first so that the threads don't conflict
  base_list = generate_base_circuits();

  //Initialize the circuit tables
  //Threading stuff
  pthread_mutex_init(&data_lock, NULL);
  pthread_mutex_init(&prnt_lock, NULL);
  pthread_cond_init(&data_ready, NULL);
  pthread_cond_init(&thrd_ready, NULL);

  pthread_t * thrds = new pthread_t[config::num_threads];
  for (i = 0; i < config::num_threads; i++) {
    pthread_create(thrds + i, NULL, &(worker_thrd), (void *)(&U));
  }
  //---------------------

  load_sequences(0, base_list, circ_table);
  pthread_mutex_lock(&data_lock);
  for (i = 1; i < config::max_seq; i++) {
    if (i % 2 == 1) {
      load_sequences(i, cliff_list, circ_table);
    } else {
      data_left = true;
      load_sequences(i, base_list, circ_table);
    }

    cout << "Looking for circuits...\n";
    clock_gettime(CLOCK_MONOTONIC, &start);

    if (i <= 1) {
      // Search the cliffords
      data_map = mp + i;
      data_left = false;
      data_circ = Circuit(1);
      data_circ.to_Rmatrix(data_mat);
      data_k = 0;
      data_avail = true;
      pthread_cond_signal(&data_ready);
      pthread_cond_wait(&thrd_ready, &data_lock);
      delete_circuit(data_circ);
    } else if (i == 2) {
      // mp[i-1] = Clifford group, so search for circuits mp[i-1]mp[i]
      data_map = mp + i;
      data_left = true;
      for (it = mp[i-1].begin(); it != mp[i-1].end(); it++) {
        data_circ = it->second;
        (it->second).to_Rmatrix(data_mat);
        for (k = 0; k < 2*pe; k += in) {
          // Write data
          data_k = k;
          data_avail = true;
          // Signal workers that data is ready
          pthread_cond_signal(&data_ready);
          pthread_cond_wait(&thrd_ready, &data_lock);
        }
      }
    } else if (i % 2 == 1) {
      // mp[i-1] = TCTC..., search for mp[i-1]mp[i-1] and mp[i]mp[i-1] 
      //  to get t-depth i-1
      data_left = false;
      for (j = 0; j < 2; j++) {
        data_map = mp + i + j - 1;
        for (it = mp[i-1].begin(); it != mp[i-1].end(); it++) {
          data_circ = it->second;
          (it->second).to_Rmatrix(data_mat);
          for (k = 0; k < 2*pe; k += in) {
            // Write data
            data_k = k;
            data_avail = true;
            // Signal workers that data is ready
            pthread_cond_signal(&data_ready);
            pthread_cond_wait(&thrd_ready, &data_lock);
          }
        }
      }
    } else {
      // mp[i-1] = CTCTC..., search for mp[i-2]mp[i] and mp[i-1]mp[i] 
      //  to get t-depth i-1
      data_map = mp + i;
      data_left = true;
      for (j = 0; j < 2; j++) {
        for (it = mp[i + j - 2].begin(); it != mp[i + j - 2].end(); it++) {
          data_circ = it->second;
          (it->second).to_Rmatrix(data_mat);
          for (k = 0; k < 2*pe; k += in) {
            // Write data
            data_k = k;
            data_avail = true;
            // Signal workers that data is ready
            pthread_cond_signal(&data_ready);
            pthread_cond_wait(&thrd_ready, &data_lock);
          }
        }
      }
    }

    /* Wait for all threads to finish */
    while (data_num > 0) {
      pthread_mutex_unlock(&data_lock);
      pthread_mutex_lock(&data_lock);
    }
    for (ti = res_list->begin(); ti != res_list->end(); ++ti) {
      (ti->second).print();
      cout << "Cost " << ti->first << "\n\n" << flush;
      delete_circuit(ti->second);
      res_list->erase(ti);
    }

    clock_gettime(CLOCK_MONOTONIC, &end);
    cout << fixed << setprecision(3);
    cout << "Time: " << (end.tv_sec + (double)end.tv_nsec/1000000000) - (start.tv_sec + (double)start.tv_nsec/1000000000) << " s\n";
    cout << "equivalent unitary vs equivalent key: " << numcorrect << " / " << numcollision << "\n";
    cout << "--------------------------------------\n" << flush;
  }

  delete [] circ_table;
  delete [] left_table;
  delete [] base_list;
  delete [] thrds;
}

void approx_search(Unitary & U) {
  int num = 0, p = 0;
  int i, j, k;
  int pe = config::mod_perms ? num_perms : 1;
  int in = config::mod_invs  ?         1 : 2;
  map_t * circ_table = new map_t[config::max_seq];
  map_t * left_table = (config::ancilla == 0) ? NULL : new map_t[config::max_seq];
  map_t * mp         = (config::ancilla == 0) ? circ_table : left_table;
  map_iter it;
  circuit_list * base_list;
  ord_circuit_list * res_list = data_res = new ord_circuit_list;
  ord_circuit_iter ti;
  struct timespec start, end;
	data_mat = Rmatrix(dim, dim);

  // Do this first so that the threads don't conflict
  base_list = generate_base_circuits();

  //Initialize the circuit tables
  //Threading stuff
  pthread_mutex_init(&data_lock, NULL);
  pthread_mutex_init(&prnt_lock, NULL);
  pthread_cond_init(&data_ready, NULL);
  pthread_cond_init(&thrd_ready, NULL);

  pthread_t * thrds = new pthread_t[config::num_threads];
  for (i = 0; i < config::num_threads; i++) {
    pthread_create(thrds + i, NULL, &(worker_thrd), (void *)(&U));
  }
  //---------------------

  load_sequences(0, base_list, circ_table);
  pthread_mutex_lock(&data_lock);
  for (i = 1; i < config::max_seq; i++) {
		if (config::ancilla == 0) {
    	load_sequences(i, base_list, circ_table);
		} else {
			load_proj(i, base_list, circ_table, left_table);
		}

    // Meet in the middle - Sequences of length 2i + {0, 1}
    data_map = mp + i;
    for (j = max(i-1, 1); j <= i; j++) {
    	cout << "Looking for circuits with depth " << 2*i - (i - j) << "...\n";
			cout << "|";
			num = 0;
			p = 0;
    	clock_gettime(CLOCK_MONOTONIC, &start);

			// Look for circuits
      for (it = mp[j].begin(); it != mp[j].end(); it++) {
        data_circ = it->second;
        (it->second).to_Rmatrix(data_mat);
        for (k = 0; k < 2*pe; k += in) {
          // Write data
          data_k = k;
          data_avail = true;
          // Signal workers that data is ready
          pthread_cond_signal(&data_ready);
          pthread_cond_wait(&thrd_ready, &data_lock);
        }
				num++;
				if (num*37 / mp[j].size() >= p) {
					cout << "=" << flush;
					p++;
				}
      }
			cout << "|\n";

			// Wait for all threads to finish 
			while (data_num > 0) {
				pthread_mutex_unlock(&data_lock);
				pthread_mutex_lock(&data_lock);
			}
			for (ti = res_list->begin(); ti != res_list->end(); ++ti) {
				(ti->second).print();
				cout << "Cost " << ti->first << "\n\n" << flush;
				delete_circuit(ti->second);
				res_list->erase(ti);
			}

			clock_gettime(CLOCK_MONOTONIC, &end);
			cout << fixed << setprecision(3);
			cout << "Time: " << (end.tv_sec + (double)end.tv_nsec/1000000000) - (start.tv_sec + (double)start.tv_nsec/1000000000) << " s\n";
			cout << "equivalent unitary vs equivalent key: " << numcorrect << " / " << numcollision << "\n";
			cout << "--------------------------------------\n" << flush;
    }
  }

  delete [] circ_table;
  if (left_table != NULL) delete [] left_table;
  delete base_list;
  delete [] thrds;
}

void mem_test(int n) {
  map_t * circ_table = new map_t[config::max_seq];
  map_t * left_table = (config::ancilla == 0) ? NULL : new map_t[config::max_seq];
  circuit_list * base_list;
  base_list = generate_base_circuits();
  load_sequences(n, base_list, circ_table);
  getchar();
}
