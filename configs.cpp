#include "configs.h"
#include <cstring>
#include <string>
#include <cstdlib>
#include <fstream>
#include <sstream>

const char * gate_names[basis_size] = {
  "I",
  "H",
  "X",  
  "Y",
  "Z",  
  "S",
  "S*",
  "T",
  "T*",
};
const char adjoint[basis_size] = {
  I,
  H,
  X,
  Y,
  Z,
  Sd,
  S,
  Td,
  T,
};
const int cnot_cost[config::num_arch] = {4,5};
const int gate_cost[config::num_arch][basis_size] = {
  {0,0,0,0,0,0,0,10,10}, 
  {0,10,0,0,0,40,40,1000,1000}
};

int num_qubits      = 0;
int num_qubits_proj = 0;
int dim             = 0;
int dim_proj        = 0;
int num_perms       = 0;
int num_weyl        = 0;

// Defaults

namespace config {
  const char * options[NUM_OPTIONS][2] = {
    {"-key-type=PROJ", "   Use projection matrices for key values"},
    {"-key-type=ACTN", "   Store the action of the unitary as its key"},
    {"-key-dim", "         Set the dimension for key generation"},
    {"-precision", "       Set the precision for approximations"},
    {"-no-phase", "        Don't mod out phase equivalences"},
    {"-no-perms", "        Don't mod out permutation equivalences"},
    {"-no-invs", "         Don't mod out inversions"},
    {"-max-seq-length", "  Set the maximum length of sequences to compute"},
    {"-no-equiv-checks", " Remove equivalence checks for key collisions"},
    {"-use-hash-map", "    Use hash table to store circuits"},
    {"-compute-tensors", " Compute rather than look up tensor products"},
    {"-hash-ring", "       Use a hash table for multiplication over R"},
    {"-tdepth", "          Search by T-depth"},
    {"-no-serialize", "    Don't store/load to/from database files"},
    {"-arch-steane", "     Use the Steane code architecture"},
    {"-arch-surface", "    Use the surface code architecture"},
    {"-approximate", "     Approximate the given unitary"},
    {"-threads", "         Set the number of threads to run on"},
    {"-space-saver", "     Turn space saving mode on"},
    {"-ancilla", "         Set the number of ancilla"},
    {"-paulis", "          Include the Pauli group in the instruction set"},
    {"-memtest", "         Run a memory test"},
		{"-frobenius", "       Turn off frobenius norm"},
    {"-rotation", "        Approximate a rotation matrix diag(1, e^i(PI * n / d))"},
    {"-help", ""}
  };
  const char circuit_file[] = "searches";

  int   key_dimension  = 1;
  key_t key_type       = PROJECTION;
  int   precision      = 1;
  bool  mod_phase      = true;
  bool  mod_perms      = true;
  bool  mod_invs       = true;
  int   max_seq        = 50;
  int   max_cliff      = 50;
  bool  check_equiv    = true;
  bool  ordered_map    = true;
  bool  tensors        = false;
  bool  hash_ring      = false;
  bool  tdepth         = false;
  bool  serialize      = true;
  Arch  architecture   = SURFACE;
  bool  approximate    = false;
  int   num_threads    = 4;
  bool  save_space     = false;
  int   ancilla        = 0;
  bool  paulis         = false;
	bool  frob_norm      = true;

  void output_config(ofstream & out) {
    out << key_dimension;
    out << (int)key_type;
    out << mod_phase;
    out << mod_perms;
    out << mod_invs;
    out << tdepth;
  }

  void input_config(ifstream & in) {
    int kd, kt;
    bool mp, mpe, mi, ms, td;
    in >> kd;
    in >> kt;
    in >> mp;
    in >> mpe;
    in >> mi;
    in >> td;
    if (kd != key_dimension || kt != (int)key_type
        || mp != mod_phase || ms != mod_perms 
        || mi != mod_invs || td != tdepth) {
      cout << "ERROR: Databases file has incompatible configuration";
      exit(1);
    }
  }
}

void init_configs(int n) {
  num_qubits = n+config::ancilla;
  num_qubits_proj = n;
  dim = 1 << num_qubits;
  dim_proj = 1 << n;
  num_perms = 1;
  for (int i = num_qubits; i > 1; i--) {
    num_perms *= i;
  }
  num_weyl = dim*dim;
}

string gen_filename(int q, int p, int d) {
  stringstream ret;
  ret << "./libraries/data" << q << "q" << p << "p" << d << "d";
  if (config::paulis) ret << "-pauli";
  return ret.str();
}

string gen_cliff_filename(int q, int d) {
  stringstream ret;
  ret << "./libraries/clifford" << q << "q" << d << "d";
  return ret.str();
}
