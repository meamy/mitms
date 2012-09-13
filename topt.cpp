#include "circuit.h"

#include <list>
#include <boost/dynamic_bitset.hpp>


// Parse a {CNOT, T} circuit into a set of phase multiples and final xor
pair<list<pair<char, boost::dynamic_bitset<> > >, boost::dynamic_bitset<> * > parse_circuit(int n, Circuit & input) {
	// Start each wire out with only its input in the sum
	boost::dynamic_bitset<> * wires = new boost::dynamic_bitset<> [n];
	for(int i = 0; i < n; i++) wires[i] = boost::dynamic_bitset<> (n, 1 << i);

	list<pair<char, boost::dynamic_bitset<> > > phases;
	list<pair<char, boost::dynamic_bitset<> > >::iterator it;
	Circuit tmp = input;
	gate cur;
	bool flg;

	while(!tmp.empty()) {
		cur = tmp.first();
		for(int i = 0; i < n; i++) {
			if (IS_C(cur[i])) {
				wires[GET_TARGET(cur[i])] ^= wires[i];
			} else if (cur[i] == T) {
				flg = false;
				for (it = phases.begin(); it != phases.end(); it++) {
					if (it->second == wires[i]) {
						it->first = (it->first + 1) % 8;
						flg = true;
					}
				}
				if (!flg) {
					phases.push_back(make_pair(1, boost::dynamic_bitset<>(wires[i])));
				}
			} else if (cur[i] == Td) {
				flg = false;
				for (it = phases.begin(); it != phases.end(); it++) {
					if (it->second == wires[i]) {
						it->first = ((it->first - 1) % 8 + 8) % 8;
						flg = true;
					}
				}
				if (!flg) {
					phases.push_back(make_pair(7, boost::dynamic_bitset<>(wires[i])));
				}
			} else if (cur[i] != I && cur[i] != X) {
				cout << "ERROR: not a {CNOT, T} circuit\n";
				phases.clear();
				delete[] wires;
				return make_pair(phases, wires);
			}
		}
		tmp = tmp.next();
	}
	
	return make_pair(phases, wires);
}

// Make triangular to determine the rank
bool is_linearly_ind(int n, const boost::dynamic_bitset<> * bits) {
	int k;

	// Make a copy of the bitset
	boost::dynamic_bitset<> * tmp = new boost::dynamic_bitset<>[n];
	for(int i = 0; i < n; i++) {
		tmp[i] = bits[i];
	}

	// Make triangular
	for(int i = 0; i < n; i++) {
		bool flg = false;
		for (int j = i; j < n; j++) {
			if (bits[j].test(i)) {
				// If we haven't yet seen a vector with bit i set...
				if (!flg) {
					// If it wasn't the first vector we tried, swap to the front
					if (j != i) swap(tmp[i], tmp[j]);
					flg = true;
				} else {
					tmp[j] ^= tmp[i];
				}
			}
		}
		if (!flg) return false;
	}
	return true;
}

// Do gaussian elimination and keep track of the xors
list<pair<int, int> > gaussian(int n, const boost::dynamic_bitset<> * bits)  {
	int i, j, k;
	list<pair<int, int> > ret;

	// Make a copy of the bitset
	boost::dynamic_bitset<> * tmp = new boost::dynamic_bitset<>[n];
	for (i = 0; i < n; i++) {
		tmp[i] = bits[i];
	}

	// Make triangular
	for (i = 0; i < n; i++) {
		bool flg = false;
		for (j = i; j < n; j++) {
			if (bits[j].test(i)) {
				// If we haven't yet seen a vector with bit i set...
				if (!flg) {
					// If it wasn't the first vector we tried, swap to the front
					if (j != i) {
						swap(tmp[i], tmp[j]);
						ret.push_back(make_pair(i, j));
						ret.push_back(make_pair(j, i));
						ret.push_back(make_pair(i, j));
					}

					flg = true;
				} else {
					tmp[j] ^= tmp[i];
					ret.push_back(make_pair(i, j));
				}
			}
		}
		if (!flg) {
			cout << "ERROR: not linearly independent\n";
			exit(1);
		}
	}

	//Finish the job
	for (i = n-1; i > 0; i--) {
		for (j = i - 1; j >= 0; j--) {
			if (bits[j].test(i)) {
				tmp[j] ^= tmp[i];
				ret.push_back(make_pair(i, j));
			}
		}
	}

	return ret;
}

// Make the input matrix a permutation and report the xors used
pair<boost::dynamic_bitset<> *, boost::dynamic_bitset<> *> 
make_per(int n, boost::dynamic_bitset<> * bits) {
	int i, j, k;
	bool flg;
	boost::dynamic_bitset<> * per = new boost::dynamic_bitset<>[n];
	boost::dynamic_bitset<> * xs  = new boost::dynamic_bitset<>[n];

	// Make a copy of the bitset
	for(i = 0; i < n; i++) {
		per[i] = bits[i];
	}

	// Make permutation
	for(i = 0; i < n; i++) {
		k = -1;

		// Set the vector with bit i
		for (j = 0; j < n; j++) {
			if (per[j].test(i)) {
				// Make sure none of the earlier bits are set
				flg = false;
				for (k = 0; k < j; k++) {
					flg = flg || per[j].test(k);
				}
				if (!flg) {
					k = j;
				}
			}
		}
		
		// If we didn't find a vector with bit i set, fail
		if (k == -1) {
			cout << "ERROR: input vectors not linearly independent\n";
			exit(1);
		}

		// Clear out all other vectors
		for (j = 0; j < n; j++) {
			if (per[j].test(i)) {
				per[j] ^= per[k];
				xs[j].flip(k);
			}
		}
	}
}

bool divide_into_bases(int n, list<pair<char, boost::dynamic_bitset<> > > lst,
	 int * m, char ** phases, boost::dynamic_bitset<> ** basis) {

	bool flg;
	int x = lst.size();
	boost::dynamic_bitset<> * tmp = new boost::dynamic_bitset<> [lst.size()];
	list<pair<char, boost::dynamic_bitset<> > >::iterator it;

	// Set up the array for permutations
	int perms[x], i;
	for (i = 0; i < x; i++) {
		perms[i] = i;
	}

	do {
		// Permute the bits according to the current permutation
		for (it = lst.begin(), i = 0; it != lst.end(); it++, i++) {
			tmp[perms[i]]  = it->second;
		}

		// Check whether they're all linearly independent
		flg = true;
		for (int i = 0; i < x / n; i++) {
			flg = flg && is_linearly_ind(n, tmp + i*n);
		}

		// If they are, quit. Otherwise, choose the next permutation
		if (flg) {
			*m = (lst.size() + n - 1) / n;
			*phases = new char [m*n];
			*basis = new boost::dynamic_bitset<> [m*n];
			for (it = lst.begin(), i = 0; it != lst.end(); it++, i++) {
				phases[perms[i]] = it->first;
				basis[perms[i]]  = it->second;
			}
			return true;
		}
	} while (next_permutation(perms, perms + x));
	return false;
}

Circuit CNOTs_to_circuit(int n, list<pair<int, int> > lst) {
	Circuit ret(lst.size());
	Circuit acc = ret;
	list<pair<int, int> >::iterator it;

	for (it = lst.begin(); it != lst.end(); it++) {
		acc[0][it->first]  = C(it->second);
		acc[0][it->second] = X;
		acc = acc.next();
	}
	return ret;
}

Circuit compute_CNOT_network(int n, const char * phases, const boost::dynamic_bitset<> * comp) {
	int i, j;
	// Compute steps to reduce to identity
	list<pair<int, int> > last  = gaussian(n, comp);

	// Do the inverse
	last.reverse();
	// Compute cnot network
	Circuit cnots = CNOTs_to_circuit(n, last);

	// Compute maximum number of T-stages
	int m = 0;
	for (i = 0; i < n; i++) m = max(m, (int)min((int)phases[i], 8 - (int)phases[i]));
	Circuit tees = Circuit(m);

	// Set the T gates
	for (i = 0; i < n; i++) {
		if (phases[i] <= 4) {
			for (j = 0; j < phases[i]; j++) tees[j][i] = T;
		} else {
			for (j = 0; j < 8 - phases[i]; j++) tees[j][i] = Td;
		}
	}

	Circuit ret = cnots.append(tees);
	delete_circuit(cnots);
	delete_circuit(tees);

	return ret;
}

Circuit compute_CNOT_network(int n, const char * phases, 
		const boost::dynamic_bitset<> * comp, const boost::dynamic_bitset<> * inp) {
	int i, j;
	list<pair<int, int> > first = gaussian(n, inp);
	// Compute steps to reduce to identity
	list<pair<int, int> > last  = gaussian(n, comp);

	// Do the inverse
	last.reverse();
	first.splice(first.end(), last);
	// Compute cnot network
	Circuit cnots = CNOTs_to_circuit(n, first);

	// Compute maximum number of T-stages
	int m = 0;
	for (i = 0; i < n; i++) m = max(m, (int)min((int)phases[i], 8 - (int)phases[i]));
	Circuit tees = Circuit(m);

	// Set the T gates
	for (i = 0; i < n; i++) {
		if (phases[i] <= 4) {
			for (j = 0; j < phases[i]; j++) tees[j][i] = T;
		} else {
			for (j = 0; j < 8 - phases[i]; j++) tees[j][i] = Td;
		}
	}

	Circuit ret = cnots.append(tees);
	delete_circuit(cnots);
	delete_circuit(tees);

	return ret;
}

Circuit parallelize(int n, int m, const char * phases, const boost::dynamic_bitset<> * bits) {
	Circuit acc(0), net, tmp;

	for (int i = 0; i < m / n; i++) {
		net = i == 0 ? compute_CNOT_network(n, phases, bits) 
								 : compute_CNOT_network(n, phases + i*n, bits + i*n, bits + (i-1)*n);
		tmp = acc;
		acc = tmp.append(net);
		delete_circuit(net);
		delete_circuit(tmp);
	}

	return acc;
}

int main() {
	num_qubits = 3;
	int n = 3;
	int m;
	Circuit x(8);
	x.circuit[0] = Td;
	x.circuit[1] = Td;
	x.circuit[2] = I;
	x.circuit[3] = X;
	x.circuit[4] = I;
	x.circuit[5] = C(0);
	x.circuit[6] = T;
	x.circuit[7] = C(2);
	x.circuit[8] = X;
	x.circuit[9] = X;
	x.circuit[10] = C(0);
	x.circuit[11] = T;
	x.circuit[12] = Td;
	x.circuit[13] = C(2);
	x.circuit[14] = X;
	x.circuit[15] = X;
	x.circuit[16] = I;
	x.circuit[17] = C(0);
	x.circuit[18] = T;
	x.circuit[19] = I;
	x.circuit[20] = Td;
	x.circuit[21] = X;
	x.circuit[22] = C(0);
	x.circuit[23] = I;
	pair<list<pair<char, boost::dynamic_bitset<> > >, boost::dynamic_bitset<> * > ret = parse_circuit(n, x);

	list<pair<char, boost::dynamic_bitset<> > >::iterator it;
	cout << "w^(";
	for (it = ret.first.begin(); it != ret.first.end(); it++) {
		if (it != ret.first.begin() && (it->first <= 4)) {
			cout << "+" << (int)(it->first);
		} else if (it->first) {
			cout << (int)(it->first - 8);
		}
		for(int i = 0; i < n; i++) {
			if (it->second[i]) cout << (char)('a' + i);
		}
	}
	cout << ")\n";

	char * phases;
	boost::dynamic_bitset<> * basis;

	divide_into_bases(n, ret.first, &m, &phases, &basis);

	Circuit tmp = parallelize(n, m, phases, basis);
	tmp.print();

	return 1;
}
