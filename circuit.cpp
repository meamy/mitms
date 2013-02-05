#include "circuit.h"
#include <cstring>
#include <iomanip>

Circuit::Circuit() { depth = 0; circuit = NULL; }
Circuit::Circuit(int d) { 
  depth = d; 
  circuit = new char[num_qubits*d]; 
  for (int i = 0; i < num_qubits*d; i++) circuit[i] = I; 
}
Circuit::Circuit(const Circuit & C) { depth = C.depth; circuit = C.circuit; }

void delete_circuit(Circuit & C) {
  delete [] C.circuit;
}

Circuit & Circuit::operator=(const Circuit & C) {
  depth = C.depth;
  circuit = C.circuit;
  return *this;
}

gate Circuit::operator[](int i) {
	return circuit + num_qubits*i;
}

Circuit & Circuit::operator++() {
  int i = 0;
  next_gate(circuit + i);
  while(is_eye(circuit+i) && (i += num_qubits) < num_qubits*depth) {
    next_gate(circuit + i);
  }
  return *this;
}

Circuit & Circuit::cliffpp() {
  int i = 0;
  next_clifford(circuit + i);
  while(is_eye(circuit+i) && (i += num_qubits) < num_qubits*depth) {
    next_clifford(circuit + i);
  }
  return *this;
}

gate Circuit::first() const {
  return circuit;
}
Circuit Circuit::next() const {
  Circuit ret;
  ret.depth = depth - 1;
  ret.circuit = circuit + num_qubits;
  return ret;
}
gate Circuit::last() const {
  return circuit + (depth-1)*num_qubits;
}

bool Circuit::empty() const {
  return (depth == 0);
}
bool Circuit::eye() const {
  for (int i = 0; i < depth*num_qubits; i++) {
    if (circuit[i] != I) return false;
  }
  return true;
}

Circuit Circuit::copy() const { 
  Circuit ret(depth);
  memcpy(ret.circuit, circuit, num_qubits*depth);
  return ret;
}

Circuit Circuit::append(const Circuit & C) const {
  Circuit ret(depth + C.depth);
  memcpy(ret.circuit, circuit, depth*num_qubits);
  memcpy(ret.circuit + depth*num_qubits, C.circuit, C.depth*num_qubits);
  return ret;
}

Circuit Circuit::reverse() const {
  Circuit ret(depth);
  for (int i = 0; i < depth; i++) {
    memcpy(ret.circuit + i*num_qubits, circuit + (depth - 1 - i)*num_qubits, num_qubits);
  }
  return ret;
}

Circuit Circuit::transform(const char * perm, bool adj) const {
  Circuit ret(depth);
  int i, j;
  for (i = 0; i < depth; i++) {
    j = adj ? depth - 1 - i : i;
    gate_transform(circuit + i*num_qubits, ret.circuit + j*num_qubits, perm, adj);
  }
  return ret;
}

Circuit Circuit::transform(int i, bool adj) const {
  const char * perm = from_lexi(i);
  return this->transform(perm, adj);
}


void Circuit::print() const {
  int i, j;
  char g;

  for (i = 0; i < num_qubits; i++) {
    for (j = (depth-1)*num_qubits + i; j >= i; j -= num_qubits) {
      g = circuit[j];
      if (IS_C(g)) {
        cout << "C(" << (int)GET_TARGET(g) + 1 << ")";
      } else {
        assert(g <= basis_size);
        cout << " " << setw(3) << left << gate_names[g];
      }
    }
    cout << "\n";
  }
}

void Circuit::to_Rmatrix(Rmatrix & U, bool adj) const {
  int i, j;
  Rmatrix V(dim, dim);

  j = adj ? depth - 1 : 0;
  gate_to_Rmatrix(circuit + j, U, adj);

  for (i = 1; i < depth; i++) {
    j = adj ? depth - 1 - i : i;
    gate_to_Rmatrix(circuit + j*num_qubits, V, adj);
    U *= V;
  }
}

void Circuit::to_Unitary(Unitary & U, bool adj) const {
  Rmatrix tmp(dim, dim);
  this->to_Rmatrix(tmp, adj);
  tmp.to_Unitary(U);
}

int Circuit::steane_cost_helper(int q) const {
  if (depth <= 0) return 0;
  else if (IS_C(circuit[q])) { 
    return cnot_cost[config::architecture] + 
      max((this->next()).steane_cost_helper(q), 
          (this->next()).steane_cost_helper(GET_TARGET(circuit[q])));
  } else if (circuit[q] == X) {
    for (int i = 0; i < num_qubits; i++) {
      if (IS_C(circuit[i]) && GET_TARGET(circuit[i]) == q) {
        return cnot_cost[config::architecture] + 
          max((this->next()).steane_cost_helper(q), 
              (this->next()).steane_cost_helper(i));
      }
    }
    return gate_cost[config::architecture][circuit[q]] + 
      (this->next()).steane_cost_helper(q);
  } else {
    return gate_cost[config::architecture][circuit[q]] + 
      (this->next()).steane_cost_helper(q);
  }
}

int Circuit::cost() const {
  int ret = 0, i;

  switch (config::architecture) {
    case config::STEANE:
      ret = this->steane_cost_helper(0);
      for (i = 1; i < num_qubits; i++) {
        ret = max(ret, this->steane_cost_helper(i));
      }
      break;
    case config::SURFACE:
      for ( i = 0; i < depth*num_qubits; i++) {
        if (IS_C(circuit[i])) {
          ret += cnot_cost[config::architecture];
        } else {
          ret += gate_cost[config::architecture][circuit[i]];
        }
      }
      break;
    default:
      break;
  }
  return ret;
}

Circuit read_circuit(istream & in) {
  Circuit ret;
  char ch[10];
  char buf[50];
  int i, j, k;
  bool flg;

  ret.depth = 1;

  for (i = config::ancilla; i < num_qubits; i++) {
    for (j = 0; j < ret.depth; j++) {
      in >> ch;
      flg = false;

      /* Check for Controlled-something */
      if (ch[0] == 'C' && ch[1] == '(') {
        k = atoi(ch + 2);
        if (k <= 0 || k > num_qubits) {
          cout << "ERROR: " << k << " is not a valid control qubit\n";
          exit(1);
        } else {
          buf[j] = C(k - 1 + config::ancilla);
          flg = true;
        }
      }

      /* Check for a single qubit gate */
      for (k = 0; k < basis_size && !flg; k++) {
        if (strcmp(ch, gate_names[k]) == 0) {
          buf[j] = k;
          flg = true;
        }
      }

      /* Other stuff */
      if (!flg) {
        cout << "ERROR: unknown gate \"" << ch << "\"";
        exit(1);
      } else if (in.peek() == '\n' && j != ret.depth-1) {
        cout << "ERROR: incomplete circuit\n";
        exit(1);
      } else if (in.peek() != '\n' && j == ret.depth-1) {
        ret.depth++;
      }
    }

    if (ret.circuit == NULL) {
      ret.circuit = new char[ret.depth*num_qubits];
    }
    for (j = 0; j < ret.depth; j++) {
      ret.circuit[j*num_qubits + i] = buf[j];
    }

    in.ignore();
  }

  /* Add ancillas */
  for (i = 0; i < config::ancilla; i++) {
    for (j = 0; j < ret.depth; j++) {
      ret.circuit[j*num_qubits + i] = I;
    }
  }

  return ret.reverse();
}

void Circuit::output(ofstream & out, int d) const {
  assert(depth >= d);
  for (int i = 0; i < d; i++) {
    output_gate(circuit + (depth - i - 1)*num_qubits, out);
  }
}
void Circuit::input (ifstream & in, int d) {
  if (d != depth) {
    if (circuit != NULL) delete [] circuit;
    depth = d;
    circuit = new char[depth];
  }

  for (int i = 0; i < depth; i++) {
    input_gate(circuit + (depth - i - 1)*num_qubits, in);
  }
}

void test_circuit() {
  Circuit A(2);

  A.circuit[0] = H;
  A.circuit[1] = C(2);
  A.circuit[2] = T;
  A.circuit[3] = X;
  A.circuit[4] = Y;
  A.circuit[5] = Z;

  Circuit C = (A.next()).reverse();
  Circuit D = (A.next()).adj();
  Circuit E = (A.next()).permute(2);
  Circuit F = C.append(D);
  
  Rmatrix a(dim, dim), b(dim, dim), c(dim, dim);
  (A.next()).to_Rmatrix(a);
  D.to_Rmatrix(b);
  a.adj(c);
  assert(c == b);
  E.to_Rmatrix(b);

  for (int i = 0; i < num_perms; i++) {
    delete_circuit(E);
    E = (A.next()).permute(i);
    E.to_Rmatrix(b);
    a.permute(c, i);
    assert(c == b);
  }

  delete_circuit(A);
  delete_circuit(C);
  delete_circuit(D);
  delete_circuit(E);
  delete_circuit(F);
}
