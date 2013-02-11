#include "search.h"
#include <cstring>

using namespace config;

Circuit * circ_search = NULL;
Unitary * unit_search = NULL;

void test_all() {
  cout << "Testing ring library......";
  test_ring();
  cout << "Complete\nTesting matrix library....";
  test_rmatrix();
  cout << "Complete\nTesting gate library......";
  test_gate();
  cout << "Complete\nTesting circuit library...";
  test_circuit();
  cout << "Complete\n";
}

void bootstrap(int n) {
  init_configs(n);
  init_ring();
  init_rmatrix();
  init_gate();
  init_util();
}

void parse_options(int argc, char *argv[]) {
  int i, tmp, n, d;
  char buf[80];

  if (argc == 1) {
    cout << "QCopt -- A tool for optimally decomposing unitaries over FT gate sets\n"
      << "Written by Matthew Amy\n"
      << "Run with QCopt [options] gate-label\n\n";
    for (int j = 0; j < NUM_OPTIONS; j++) {
      cout << options[j][0] << "  " << options[j][1] << "\n";
    }
    exit(0);
  }

  for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			if (strcmp(argv[i], options[0][0]) == 0) {
				key_type = PROJECTION;
			} else if (strcmp(argv[i], options[1][0]) == 0) {
				key_type = ACTION;
			} else if (strcmp(argv[i], options[2][0]) == 0) {
				if (i >= argc || (tmp = atoi(argv[i+1])) < 1) {
					cout << "Key dimension must be greater than 0\n";
					exit(1);
				} else {
					key_dimension = tmp;
					i++;
				}
			} else if (strcmp(argv[i], options[3][0]) == 0) {
				if (i >= argc || (tmp = atoi(argv[i+1])) < 1) {
					cout << "Precision must be greater than 0\n";
					exit(1);
				} else {
					precision = pow(10, -tmp);
					i++;
				}
			} else if (strcmp(argv[i], options[4][0]) == 0) {
				mod_phase = false;
			} else if (strcmp(argv[i], options[5][0]) == 0) {
				mod_perms = false;
			} else if (strcmp(argv[i], options[6][0]) == 0) {
				mod_invs = false;
			} else if (strcmp(argv[i], options[7][0]) == 0) {
				if (i >= argc || (tmp = atoi(argv[i+1])) < 1) {
					cout << "Maximum sequence length must be at least 1\n";
					exit(1);
				} else {
					max_seq = tmp;
					i++;
				}
			} else if (strcmp(argv[i], options[8][0]) == 0) {
				check_equiv = false;
			} else if (strcmp(argv[i], options[9][0]) == 0) {
				ordered_map = false;
			} else if (strcmp(argv[i], options[10][0]) == 0) {
				tensors = true;
			} else if (strcmp(argv[i], options[11][0]) == 0) {
				hash_ring = true;
			} else if (strcmp(argv[i], options[12][0]) == 0) {
				tdepth = true;
				mod_invs = false;
				serialize = false;
			} else if (strcmp(argv[i], options[13][0]) == 0) {
				serialize = false;
			} else if (strcmp(argv[i], options[14][0]) == 0) {
				architecture = STEANE;
			} else if (strcmp(argv[i], options[15][0]) == 0) {
				architecture = SURFACE;
			} else if (strcmp(argv[i], options[16][0]) == 0) {
				approximate = true;
			} else if (strcmp(argv[i], options[17][0]) == 0) {
				if (i >= argc || (tmp = atoi(argv[i+1])) < 1) {
					cout << "At least 1 thread required to run 0\n";
					exit(1);
				} else {
					num_threads = tmp;
					i++;
				}
			} else if (strcmp(argv[i], options[18][0]) == 0) {
				save_space = true;
			} else if (strcmp(argv[i], options[19][0]) == 0) {
				if (i >= argc || (tmp = atoi(argv[i+1])) < 0) {
					cout << "Specify a non-negative number of ancillary qubits\n";
					exit(1);
				} else {
					ancilla = tmp;
					i++;
				}
			} else if (strcmp(argv[i], options[20][0]) == 0) {
				paulis = true;
			} else if (strcmp(argv[i], options[21][0]) == 0) {
				if (i >= argc || (tmp = atoi(argv[i+1])) < 0) {
					cout << "Specify a non-negative number of qubits\n";
					exit(1);
				} else {
					bootstrap(tmp);
				}
				if (i >= argc-1 || (tmp = atoi(argv[i+2])) < 0) {
					cout << "Specify a non-negative depth\n";
					exit(1);
				} else {
					mem_test(tmp);
					exit(0);
				}
			} else if (strcmp(argv[i], options[22][0]) == 0) {
				frob_norm = false;
			} else if (strcmp(argv[i], options[23][0]) == 0) {
				if (i >= argc) {
					cout << "Specify an integer numerator\n";
					exit(1);
				} else if (i >= argc-1 || (d = atoi(argv[i+2])) < 0) {
					cout << "Specify a non-zero integer denominator\n";
					exit(1);
				} else {
          n = atoi(argv[i+1]);
          bootstrap(1);
          unit_search = new Unitary(dim, dim);
          (*unit_search)(0, 0) = LaComplex(1, 0);
          (*unit_search)(0, 1) = LaComplex(0, 0);
          (*unit_search)(1, 0) = LaComplex(0, 0);
          (*unit_search)(1, 1) = LaComplex(polar(1.0, PI * n / d));
          i += 2;
				}
			} else if (strcmp(argv[i], options[24][0]) == 0) {
				cout << "QCopt -- A tool for optimally decomposing unitaries over FT gate sets\n"
					<< "Written by Matthew Amy\n"
					<< "Run with QCopt [options] gate-label\n\n";
				for (int j = 0; j < NUM_OPTIONS; j++) {
					cout << options[j][0] << "  " << options[j][1] << "\n";
				}
				exit(0);
			} else {
				cout << "Unrecognized option \"" << argv[i] << "\"\n";
			}
		} else {
			ifstream in;
			in.open(circuit_file, iostream::in);
      bool flag = false;
			while(!in.eof()) {
				in >> buf;
				if (strcmp(argv[i], buf) == 0) {
          flag = true;
					in >> tmp;
					if (tmp <= 0) {
						cout << "ERROR: invalid circuit description\n";
						exit(1);
					}
					while(isspace(in.peek())) {
						in.ignore();
					}

					bootstrap(tmp);
          circ_search  = new Circuit;
					*circ_search = read_circuit(in);
				}
			}
			if (!flag) cout << "No circuit \"" << argv[i] << "\" found\n";
		}
  }
}

int main(int argc, char * argv[]) {
  cout << "\n";
  parse_options(argc, argv);
  if (circ_search != NULL) {
    Rmatrix U(dim, dim), V(dim, dim);
    circ_search->to_Rmatrix(V);
    V.submatrix(0, 0, dim_proj, dim_proj, U);
    cout << "Searching for U = \n";
    U.print();
    if (!approximate) {
      if (tdepth) exact_search_tdepth(U);
      else exact_search(U);
    } else {
			approx_search(U);
		}
  } else if (unit_search != NULL) {
    cout << "Searching for U = \n" << *unit_search;
    approx_search(*unit_search);
  }
  cout << "\n";
  
  return 0;
}
