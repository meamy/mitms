#include "search.h"
#include <string>
#include <iomanip>
#include <pthread.h>
#include "vptree.h"

/* ------------- For use with approximate searching, for some reason I don't remember */
class map_value_iter : public std::iterator<bidirectional_iterator_tag, Circuit> {
	private:
		map_iter it;
	public:
		map_value_iter(map_iter x) : it(x) { }
		map_value_iter(const map_value_iter & mit) : it(mit.it) { }
		map_value_iter & operator++()    { ++it; return *this; }
		map_value_iter   operator++(int) { map_value_iter tmp(*this); operator++(); return tmp; }
		bool operator==(const map_value_iter & rhs) { return it == rhs.it; }
		bool operator!=(const map_value_iter & rhs) { return it != rhs.it; }
		Circuit & operator*() { return it->second; }
};

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

/* -------------------------------- */

bool check_it(const Circuit & x, const Circuit & y, const Rmatrix & target) {
  Rmatrix tmp1(dim, dim), tmp2(dim, dim);
  x.to_Rmatrix(tmp1);
  y.to_Rmatrix(tmp2);

  tmp1 *= tmp2;
  if (dim == dim_proj) {
    return target.phase_eq(tmp1);
  } else {
    tmp1.submatrix(0, 0, dim_proj, dim_proj, tmp2);
    return target.phase_eq(tmp2);
  }

}

void * worker_thrd(void * arg) {
  Circuit circ, tmp_circ, tmp_circ2, tmp_circ3, ans;
  int k, i, j, cst;
  bool left_multiply;
  map_t * mp;
  Rmatrix V(dim, dim), W(dim_proj, dim), U = *((Rmatrix *)arg);
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

  load_sequences(0, base_list, circ_table, NULL);
  pthread_mutex_lock(&data_lock);
  for (i = 1; i < config::max_seq; i++) {
		if (config::ancilla == 0) {
    	load_sequences(i, base_list, circ_table, NULL);
		} else {
			load_sequences(i, base_list, left_table, circ_table);
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
        int xxx = num*37 / (mp[j].size());
				if (xxx >= p) {
          for (int yyy = p; yyy <= ceil(xxx); yyy++) {
					  cout << "=" << flush;
					  p++;
          }
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

  load_sequences(0, base_list, circ_table, NULL);
  pthread_mutex_lock(&data_lock);
  for (i = 1; i < config::max_seq; i++) {
    if (i % 2 == 1) {
      load_sequences(i, cliff_list, circ_table, NULL);
    } else {
      data_left = true;
      load_sequences(i, base_list, circ_table, NULL);
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
    cout << "--------------------------------------\n" << flush;
  }

  delete [] circ_table;
  delete [] left_table;
  delete [] base_list;
  delete [] thrds;
}

//------------------------------------------------Approximations

double circuit_dist(const Circuit & a, const Circuit & b) {
	Rmatrix U(dim, dim), V(dim, dim);
	a.to_Rmatrix(U);
	b.to_Rmatrix(V);
	return dist(U, V);
}

// Create a projector onto [0,1] that closes over the circuit's unitary
class circuit_closure {
	private:
		Rmatrix p;
	public:
		circuit_closure(const Circuit & circ) : p(dim, dim) { circ.to_Rmatrix(p); }
		circuit_closure(const Rmatrix & mat) : p(mat) { }
		double operator()(const Circuit & circ) const { 
			Rmatrix q(dim, dim);
			circ.to_Rmatrix(q);
			return dist(p, q);
		}
};
typedef VPTree<Circuit, circuit_dist, circuit_closure> NNtree;
typedef NNtree::iterator             NN_iter;
typedef NNtree::const_iterator const_NN_iter;

void approx_search(Rmatrix & U) {
  int num = 0, p = 0;
  int i, j, k;
  int pe = config::mod_perms ? num_perms : 1;
  int in = config::mod_invs  ?         1 : 2;
  map_t  * circ_table = new map_t [config::max_seq];
	NNtree * NN_table   = new NNtree[config::max_seq];
  NN_iter it;
  circuit_list * base_list;
  struct timespec start, end;
	Rmatrix mat(dim, dim), V(dim, dim);
	const Circuit * ans;
	Circuit tmp_circ, tmp_circ2, tmp_circ3;
	double epsilon = config::precision;
  Canon * canon_form;
  struct triple * trip;

  // Do this first so that the threads don't conflict
  base_list = generate_base_circuits();

  load_sequences(0, base_list, circ_table, NULL);
  for (i = 1; i < config::max_seq; i++) {
   	load_sequences(i, base_list, circ_table, NULL);
		NN_table[i].build_tree(
				map_value_iter(circ_table[i].begin()),
			 	map_value_iter(circ_table[i].end()),
			 	circ_table[i].size()
		);

    // Meet in the middle - Sequences of length 2i + {0, 1}
    for (j = max(i-1, 1); j <= i; j++) {
    	cout << "Looking for circuits with depth " << 2*i - (i - j) << "...\n";
			cout << "|";
			num = 0;
			p = 0;
    	clock_gettime(CLOCK_MONOTONIC, &start);

			// Look for circuits
      for (it = NN_table[j].begin(); it != NN_table[j].end(); it++) {
        it->to_Rmatrix(mat);
        for (k = 0; k < 2*pe; k += in) {

					if (k % 2 == 0) {
						mat.permute_adj(V, k/2);
					} else {
						mat.permute(V, k/2);
					}

					// Perform the search
					canon_form = canonicalize(U*V);

					trip = &(canon_form->front());
					ans = NN_table[i].nearest_neighbour(circuit_closure(trip->mat), &epsilon);
					if (ans != NULL) {
						// Generate the two circuit halves
						tmp_circ = (k == 0) ? *it : it->transform(k/2, k % 2 == 1);
						tmp_circ2 = ans->transform(-(trip->permutation), trip->adjoint);
						tmp_circ3 = tmp_circ2.append(tmp_circ);

						pthread_mutex_lock(&prnt_lock);
						tmp_circ3.print();
						cout << scientific << epsilon << "\n" << flush;
						pthread_mutex_unlock(&prnt_lock);

						if (k != 0) delete_circuit(tmp_circ);
						delete_circuit(tmp_circ2);
						delete_circuit(tmp_circ3);
					}
					canon_form->clear();
					delete canon_form;

        }
				num++;
        int xxx = num*37 / (NN_table[j].size());
				if (xxx >= p) {
          for (int yyy = p; yyy <= ceil(xxx); yyy++) {
					  cout << "=" << flush;
					  p++;
          }
				}
      }
			cout << "|\n";

			clock_gettime(CLOCK_MONOTONIC, &end);
			cout << fixed << setprecision(3);
			cout << "Time: " << (end.tv_sec + (double)end.tv_nsec/1000000000) - (start.tv_sec + (double)start.tv_nsec/1000000000) << " s\n";
			cout << "--------------------------------------\n" << flush;
    }
  }

  delete [] circ_table;
  delete base_list;
}

void approx_search(Unitary & U) {
	/*
  int num = 0, p = 0;
  int i, j, k;
  int pe = config::mod_perms ? num_perms : 1;
  int in = config::mod_invs  ?         1 : 2;
  map_t  * circ_table = new map_t [config::max_seq];
	NNtree * NN_table   = new NNtree[config::max_seq];
  NN_iter it;
  circuit_list * base_list;
  struct timespec start, end;
	Rmatrix mat(dim, dim), V(dim, dim);
	Unitary tmp_V(dim, dim);
	const Circuit * ans;
	Circuit tmp_circ, tmp_circ2, tmp_circ3;
	double epsilon = config::precision;
  Canon * canon_form;
  struct triple * trip;

  // Do this first so that the threads don't conflict
  base_list = generate_base_circuits();

  load_sequences(0, base_list, circ_table);
  for (i = 1; i < config::max_seq; i++) {
   	load_sequences(i, base_list, circ_table);
		NN_table[i].build_tree(
				map_value_iter(circ_table[i].begin()),
			 	map_value_iter(circ_table[i].end()),
			 	circ_table[i].size()
		);

    // Meet in the middle - Sequences of length 2i + {0, 1}
    for (j = max(i-1, 1); j <= i; j++) {
    	cout << "Looking for circuits with depth " << 2*i - (i - j) << "...\n";
			cout << "|";
			num = 0;
			p = 0;
    	clock_gettime(CLOCK_MONOTONIC, &start);

			// Look for circuits
      for (it = NN_table[j].begin(); it != NN_table[j].end(); it++) {
        it->to_Rmatrix(mat);
        for (k = 0; k < 2*pe; k += in) {

					if (k % 2 == 0) {
						mat.permute_adj(V, k/2);
					} else {
						mat.permute(V, k/2);
					}

					// Perform the search
					V.to_Unitary(tmp_V);
					canon_form = canonicalize(U*tmp_V);

					trip = &(canon_form->front());
					ans = NN_table[i].nearest_neighbour(circuit_closure(trip->mat), &epsilon);
					if (ans != NULL) {
						// Generate the two circuit halves
						tmp_circ = (k == 0) ? *it : it->transform(k/2, k % 2 == 1);
						tmp_circ2 = ans->transform(-(trip->permutation), trip->adjoint);
						tmp_circ3 = tmp_circ2.append(tmp_circ);

						pthread_mutex_lock(&prnt_lock);
						tmp_circ3.print();
						cout << scientific << epsilon << "\n" << flush;
						pthread_mutex_unlock(&prnt_lock);

						if (k != 0) delete_circuit(tmp_circ);
						delete_circuit(tmp_circ2);
						delete_circuit(tmp_circ3);
					}
					canon_form->clear();
					delete canon_form;

        }
				num++;
				if (num*37 / NN_table[j].size() >= p) {
					cout << "=" << flush;
					p++;
				}
      }
			cout << "|\n";

			clock_gettime(CLOCK_MONOTONIC, &end);
			cout << fixed << setprecision(3);
			cout << "Time: " << (end.tv_sec + (double)end.tv_nsec/1000000000) - (start.tv_sec + (double)start.tv_nsec/1000000000) << " s\n";
			cout << "--------------------------------------\n" << flush;
    }
  }

  delete [] circ_table;
  delete base_list;
	*/
}

void mem_test(int n) {
  map_t * circ_table = new map_t[config::max_seq];
  map_t * left_table = (config::ancilla == 0) ? NULL : new map_t[config::max_seq];
  circuit_list * base_list;
  base_list = generate_base_circuits();
  load_sequences(n, base_list, circ_table, NULL);
  getchar();
}
