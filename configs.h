
/*--------------------------------------------------------------------
MITMS - Meet-in-the-middle quantum circuit synthesis
Copyright (C) 2013  Matthew Amy and The University of Waterloo,
Institute for Quantum Computing, Quantum Circuits Group

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Author: Matthew Amy
---------------------------------------------------------------------*/

/* Configuration settings and necessary definitions for program */

#ifndef CONFIGS
#define CONFIGS

#define LA_COMPLEX_SUPPORT
#define PI M_PI
#define NUM_OPTIONS 26

#include <iostream>
#include <fstream>

/* If the c++ std lib implementation does not have unordered_map.h, set this to false */
#define HAS_HASH true
#if HAS_HASH
  #include <unordered_map>
#else
  #include <map>
#endif

/* Basis state projections on one qubit */
#define PROJ(x, y)     (basis_size + 2*x + y)
#define IS_PROJ(x)  (basis_size <= x && x <= 0x3e)
#define GET_PROJ(x) (x - basis_size)

/* Controls. Highest order bit defines a control, the rest of the
   byte specifies the target */
#define C(x)          (x | 0x40)
#define IS_C(x)       (x & 0x40)
#define GET_TARGET(x) (x & 0x3F)

using namespace std;

template< class Key, class Typ, class Hash, class Pred >
struct hash_table {
#if HAS_HASH
  typedef unordered_map<Key, Typ, Hash, Pred> t;
#else
  typedef multimap<Key, Typ, Pred> t;
#endif
};

/* 1 qubit basis gates */
enum basis_gates {
  I = 0,
  H = 1,
  X = 2,
  Y = 3,
  Z = 4,
  S = 5,
  Sd= 6,
  T = 7,
  Td= 8,
  basis_size = 9
};

/* Gate constants */
extern const char * gate_names[basis_size];
extern const char   adjoint[basis_size];
extern const int    gate_cost[][basis_size];
extern const int    cnot_cost[];

/* Runtime configurations */
extern int num_qubits;
extern int num_qubits_proj;
extern int dim;
extern int dim_proj;
extern int num_perms;
extern int num_weyl;

namespace config {
  enum Arch { 
    STEANE = 0, 
    SURFACE = 1,
    num_arch = 2
  };

  enum key_t {
    PROJECTION,
    ACTION
  };

  extern const char * options[NUM_OPTIONS][2];
  extern const char circuit_file[];

  /* Command line configurations */
  extern int    key_dimension; // Number of vectors used in key construction
  extern key_t  key_type;       // Define the type of key to use
  extern double precision; // Precision
  extern bool   mod_perms; // Whether we mod out permutations
  extern bool   mod_invs; // Whether we mod out inverse
  extern bool   mod_phase;      // Whether we mod out phase
  extern int    max_seq; // Max circuit length
  extern int    max_cliff; // Max depth we try to find clifford circuits
  extern bool   check_equiv; // Whether we check to make sure two circuits are equiv
  extern bool   ordered_map; // Whether we should use an ordered map
  extern bool   tensors; // Whether to store gates as tensor products of gates
  extern bool   hash_ring; // Whether we want to use a hash table to perform lookups in the ring
  extern bool   tdepth; // Whether we want to search by T-depth
  extern bool   serialize; // Whether to serialize the generated circuits
  extern Arch   architecture; // Which fault tolerant architecture
  extern bool   approximate; // Whether we're approximating a unitary
  extern int    num_threads; // How many threads to run
  extern bool   save_space; // Whether to use some space saving techniques
  extern int    ancilla; // Number of ancilla qubits to use
  extern bool   paulis; // Whether to include the paulis in the instruction set
  extern bool   frob_norm; // Whether to use the frobenius norm

  /* Configuration tools */
  void output_config(ofstream & out);
  void input_config (ifstream & in);
}

void init_configs(int n);
string gen_filename(int q, int p, int d);
string gen_cliff_filename(int q, int d);

#endif
