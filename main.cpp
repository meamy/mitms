
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

void ignore_white(istream& in) {
  while (in.peek() == ' ' || in.peek() == ';') in.ignore();
}

// Doesn't handle ancillae
Circuit read_dotqc(istream& in) {
  int i, j, n = 0, m;
  string buf, tmp;
  list<string> names;
  list<string>::iterator nit;
  list<pair<string, list<int> > > circ;
  list<int> q_list;

  // Inputs
  while (buf != ".v") in >> buf;
  ignore_white(in);
  while(in.peek() != '\n' && in.peek() != '\r') {
    in >> buf;
    names.push_back(buf);
    ignore_white(in);
  }

  // Circuit
  while (buf != "BEGIN") in >> buf;
  in >> tmp;
  while (tmp != "END") {
    q_list.clear();
    // Build up a list of the applied qubits
    ignore_white(in);
    while (in.peek() != '\n' && in.peek() != '\r' && in.peek() != ';') {
      in >> buf;
      int pos = buf.find(';');
      if (pos != string::npos) {
        for (i = buf.length() - 1; i > pos; i--) {
          in.putback(buf[i]);
        }
        in.putback('\n');
        buf.erase(pos, buf.length() - pos);
      }
      for (nit = names.begin(), i = 0; nit != names.end(); nit++, i++) {
        if (*nit == buf) q_list.push_back(i);
      }
      ignore_white(in);
    }
    if (tmp == "TOF") tmp = "tof";
    circ.push_back(make_pair(tmp, q_list));
    in >> tmp;
  }

  bootstrap(names.size());
  Circuit ret(circ.size());
  list<pair<string, list<int> > >::iterator it;
  list<int>::iterator ti;
  for (it = circ.begin(), i = 0; it != circ.end(); it++, i++) {
    ti = it->second.end();
    ti--;
    if (it->first == "T") ret.circuit[i*num_qubits + *ti] = T;
    else if (it->first == "T*") ret.circuit[i*num_qubits + *ti] = Td;
    else if (it->first == "P") ret.circuit[i*num_qubits + *ti] = S;
    else if (it->first == "P*") ret.circuit[i*num_qubits + *ti] = Sd;
    else if (it->first == "X") ret.circuit[i*num_qubits + *ti] = X;
    else if (it->first == "Y") ret.circuit[i*num_qubits + *ti] = Y;
    else if (it->first == "Z") ret.circuit[i*num_qubits + *ti] = Z;
    else if (it->first == "tof") ret.circuit[i*num_qubits + *ti] = X;
    else if (it->first == "H") ret.circuit[i*num_qubits + *ti] = H;
    for (list<int>::iterator iti = it->second.begin(); iti != ti; iti++) {
      ret.circuit[i*num_qubits + *iti] = C(*ti + config::ancilla);
    }
  }
  return ret.reverse();
}

void parse_options(int argc, char *argv[]) {
  int i, tmp, n, d;
  char buf[80];

  if (argc == 1) {
    cout << "MITMS -- A tool for optimally decomposing unitaries over FT gate sets\n"
      << "Written by Matthew Amy\n"
      << "Run with mitms [options] <circuit-label>\n\n";
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
          precision = pow(10.0, -tmp);
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
        Circuit tmp = read_dotqc(cin);
        tmp.print();
        cout << "\n";
        Rmatrix U(dim, dim);
        tmp.to_Rmatrix(U);
        U.print();
      } else if (strcmp(argv[i], options[25][0]) == 0) {
        cout << "MITMS -- A tool for optimally decomposing unitaries over FT gate sets\n"
          << "Written by Matthew Amy\n"
          << "Run with mitms [options] <circuit-label>\n\n";
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
