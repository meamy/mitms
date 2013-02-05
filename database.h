#ifndef DATABASE
#define DATABASE

#include "util.h"

/* Circuit databases */
typedef multimap<hash_t, Circuit, cmp_hash> map_t;
typedef map_t::iterator                     map_iter;
typedef map_t::const_iterator         const_map_iter;
typedef pair    <hash_t, Circuit>           map_elt;

/* Structure for returning results of searches */
typedef struct result {
  bool first;
  Circuit second;
  map_iter third;
  result() { }
  result(bool b, Circuit c) { first = b; second = c; }
  result(bool b, Circuit c, map_iter it) { first = b; second = c; third = it; }
} result;

/* Circuit lists */
typedef list<Circuit>          circuit_list;
typedef circuit_list::iterator circuit_iter;

/* Ordered Circuit lists */
typedef multimap<int, Circuit>     ord_circuit_list;
typedef ord_circuit_list::iterator ord_circuit_iter;
typedef pair    <int, Circuit>     ord_circuit_pair;

/* List of the clifford group circuits */
extern circuit_list * cliff_list;

result find_unitary(const hash_t & key, const Unitary & U, map_t & map);
result find_unitary(const hash_t & key, const Rmatrix & U, map_t & map);
circuit_list * generate_base_circuits();
void load_sequences(int i, circuit_list * L, map_t * circ_table, map_t * unproj);

#endif
