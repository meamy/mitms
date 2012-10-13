#include "circuit.h"

#include <list>
#include <set>
#include <boost/dynamic_bitset.hpp>

// Data structures for a {CNOT, T} circuit
typedef pair<char, boost::dynamic_bitset<> > exponent;
struct cnot_t_circuit {
	int                       num_inputs;
	list<exponent> 						phase_exponents;
	boost::dynamic_bitset<> * output_functions;

	void print() {
		int i, j;
		list<exponent>::iterator it;

		cout << "U|";
		for (i = 0; i < num_inputs; i++) {
			cout << (char)('a' + i);
		}
		cout << "> --> w^(";

		// Print the phase exponents
		for (it = phase_exponents.begin(); it != phase_exponents.end(); it++) {
			if (it != phase_exponents.begin() && (it->first <= 4)) {
				cout << "+" << (int)(it->first);
			} else if (it->first) {
				cout << (int)(it->first - 8);
			}
			for(int i = 0; i < num_inputs; i++) {
				if (it->second[i]) cout << (char)('a' + i);
			}
		}
		cout << ")|";

		// Print the output functions
		for (i = 0; i < num_inputs; i++) {
			cout << "(";
			for (j = 0; j < num_inputs; j++) {
				if (output_functions[i].test(j)) {
					cout << (char)('a' + j);
				}
			}
			cout << ")";
		}
		cout << ">\n";
	}
};

// Parse a {CNOT, T} circuit
cnot_t_circuit parse_circuit(int n, const Circuit & input) {
	cnot_t_circuit ret;
	ret.num_inputs = n;
	boost::dynamic_bitset<> * wires = ret.output_functions = new boost::dynamic_bitset<> [n];

	// Start each wire out with only its input in the sum
	for(int i = 0; i < n; i++) wires[i] = boost::dynamic_bitset<> (n, 1 << i);

	list<exponent>::iterator it;
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
				for (it = ret.phase_exponents.begin(); it != ret.phase_exponents.end(); it++) {
					if (it->second == wires[i]) {
						it->first = (it->first + 1) % 8;
						flg = true;
					}
				}
				if (!flg) {
					ret.phase_exponents.push_back(make_pair(1, boost::dynamic_bitset<>(wires[i])));
				}
			} else if (cur[i] == Td) {
				flg = false;
				for (it = ret.phase_exponents.begin(); it != ret.phase_exponents.end(); it++) {
					if (it->second == wires[i]) {
						it->first = ((it->first - 1) % 8 + 8) % 8;
						flg = true;
					}
				}
				if (!flg) {
					ret.phase_exponents.push_back(make_pair(7, boost::dynamic_bitset<>(wires[i])));
				}
			} else if (cur[i] != I && cur[i] != X) {
				cout << "ERROR: not a {CNOT, T} circuit\n";
				ret.phase_exponents.clear();
				delete[] ret.output_functions;
				return ret;
			}
		}
		tmp = tmp.next();
	}
	
	return ret;
}

// Make triangular to determine the rank
int compute_rank(int m, int n, const boost::dynamic_bitset<> * bits) {
	int k;
	int ret = 0;

	// Make a copy of the bitset
	boost::dynamic_bitset<> * tmp = new boost::dynamic_bitset<>[m];
	for(int i = 0; i < m; i++) {
		tmp[i] = bits[i];
	}

	// Make triangular
	for(int i = 0; i < n; i++) {
		bool flg = false;
		for (int j = i; j < m; j++) {
			if (tmp[j].test(i)) {
				// If we haven't yet seen a vector with bit i set...
				if (!flg) {
					// If it wasn't the first vector we tried, swap to the front
					if (j != i) swap(tmp[i], tmp[j]);
					flg = true;
					ret++;
				} else {
					tmp[j] ^= tmp[i];
				}
			}
		}
	}
	delete [] tmp;
	return ret;
}

bool test_oracle(const boost::dynamic_bitset<> * bits, const set<int> lst) {
	set<int>::const_iterator it;
	int i;
	boost::dynamic_bitset<> * tmp = new boost::dynamic_bitset<>[lst.size()];
	cout << "Testing set {";
	for (i = 0, it = lst.begin(); it != lst.end(); it++, i++) {
		cout << *it << ",";
		tmp[i] = bits[*it];
	}

	cout << "}\n" << flush;
	return (compute_rank(lst.size(), 3, tmp) == lst.size());
}

// Do gaussian elimination and keep track of the xors
list<pair<int, int> > gaussian(int m, int n, const boost::dynamic_bitset<> * bits)  {
	int i, j, k;
	list<pair<int, int> > ret;

	// Make a copy of the bitset
	boost::dynamic_bitset<> * tmp = new boost::dynamic_bitset<>[m];
	for (i = 0; i < m; i++) {
		tmp[i] = bits[i];
	}

	// Make triangular
	for (i = 0; i < n; i++) {
		bool flg = false;
		for (j = i; j < m; j++) {
			if (tmp[j].test(i)) {
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
			cout << "ERROR: not full rank\n";
			exit(1);
		}
	}

	//Finish the job
	for (i = n-1; i > 0; i--) {
		for (j = i - 1; j >= 0; j--) {
			if (tmp[j].test(i)) {
				tmp[j] ^= tmp[i];
				ret.push_back(make_pair(i, j));
			}
		}
	}

	delete [] tmp;
	return ret;
}

// Given a number of qubits, target rank, an incomplete basis,
//	add vectors to the basis to make it full rank
void make_full_rank(int m, int n, boost::dynamic_bitset<> * basis) {
	int k = m - 1;

	// Make a copy of the bitset
	boost::dynamic_bitset<> * tmp = new boost::dynamic_bitset<>[m];
	for(int i = 0; i < m; i++) {
		tmp[i] = basis[i];
	}

	// Make triangular
	for(int i = 0; i < n; i++) {
		bool flg = false;
		for (int j = i; j < m; j++) {
			if (tmp[j].test(i)) {
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
		if (!flg) basis[k--].set(i);
	}

	delete [] tmp;
}

// Partition a matroid
typedef list< set< int > > partitioning;

template <class t>
partitioning partition_matroid(int n, t * elts, bool oracle(const t *, const set<int>)) {
	partitioning ret;
	partitioning::iterator Si;
	set<int>::iterator yi;

	// The node q contains a queue of paths and an iterator to each node's location.
	//	Each path's first element is the element we grow more paths from.
	//	If x->y is in the path, then we can replace x with y.
	list <list <set<int>::iterator> > node_q;
	list <set<int>::iterator>::iterator path;
	bool marked[n], flag;
	int current, tmp;
	list <set<int>::iterator> tmp_list;

	set<int> hack;
	set<int> * newset;

	// For each element of the matroid
	for (int i = 0; i < n; i++) {
		cout << "Adding " << i << " to partition: " << flush;
		for (Si = ret.begin(); Si != ret.end(); Si++) {
			cout << "{";
			for (yi = Si->begin(); yi != Si->end(); yi++) {
				cout << *yi << ",";
			}
			cout << "}";
		}
		cout << "\n" << flush;

		// Build the inverse directed graph Gx breadth-first and short circuit
		//	when we find a path to a partition
		yi = hack.insert(hack.begin(), i);
		node_q.push_back(list<set<int>::iterator>(1,yi));
		marked[i] = true;
		flag = false;

		// BFS loop
		while (!node_q.empty() && !flag) {
			// The matroid element (as an index to elts) that we're currently considering
			current = *((node_q.front()).front());
			cout << "Current node: " << current << "\n" << flush;

			// For each independent set, check if we can either add CURRENT in which case 
			//	we do so and make all the changes. Otherwise, we check if any elements
			//	can be removed to make it independent, and if so add that element as a new
			//	node in the queue
			for (Si = ret.begin(); Si != ret.end() && !flag && Si->find(current) == Si->end(); Si++) {

				// Add CURRENT to Si. If Si is independent, leave it, otherwise we'll have to remove it
				Si->insert(current);
				if (oracle(elts, *Si)) {
					cout << "Path found: Si";
					// We have the shortest path to a partition, so make the changes:
					//	For each x->y in the path, remove x from its partition and add y
					for (path = (node_q.front()).begin(); path != --((node_q.front()).end()); ) {
						cout << "-->" << **(path);
						*(path) = *(++path);
					}
					cout << "-->" << current << "\n" << flush;
					flag = true;
				} else {
					// For each element of Si, if removing it makes an independent set, add it to the queue
					for (yi = Si->begin(); yi != Si->end(); yi++) {
						// Only consider adding yi if it's not already in the graph
						if (!marked[*yi]) {
							// Take yi out
							tmp = *yi;
							yi = Si->erase(yi);
							if (oracle(elts, *Si)) {
								cout << "Adding edge " << tmp << "-->" << current << "\n" << flush;
								// Put yi back in
								yi = Si->insert(yi, tmp);
								marked[*yi] = true;
								// Add yi to the end of CURRENT's path and insert
								tmp_list = node_q.front(); // Copy the path
								tmp_list.push_front(yi);
								node_q.push_back(tmp_list);
							} else yi = Si->insert(yi, tmp);
						}
					}
					// Remove CURRENT from Si
					Si->erase(current);
				}
			}

			node_q.pop_front();
		}

		// We were unsuccessful trying to edit the current partitions
		if (!flag) {
			newset = new set<int>;
			newset->insert(i);
			ret.push_front(*newset);
		}

		// Reset everything
		node_q.clear();
		for (int j = 0; j <= i; j++) {
			marked[j] = false;
		}

	}
	return ret;
}

// Given a number of qubits, target rank, and a list of phase/vector pairs,
//  optimally divide the phase/vector pairs into k collections
bool divide_into_bases(int m, int n, list<pair<char, boost::dynamic_bitset<> > > lst,
	 int * k, char ** phases, boost::dynamic_bitset<> ** basis) {

	bool flg;
	int num = lst.size();
	boost::dynamic_bitset<> * tmp = new boost::dynamic_bitset<> [num];
	list<pair<char, boost::dynamic_bitset<> > >::iterator it;

	// Set up the array for permutations
	int perms[num], i;
	for (i = 0; i < num; i++) {
		perms[i] = i;
	}

	do {
		// Permute the bits according to the current permutation
		for (it = lst.begin(), i = 0; it != lst.end(); it++, i++) {
			tmp[perms[i]]  = it->second;
		}

		// Check whether they're all full rank, possibly except for the last one
		flg = true;
		for (int i = 0; i < num / m; i++) {
			flg = flg && (compute_rank(m, n, tmp + i*m) == n);
		}

		// If they are, quit. Otherwise, choose the next permutation
	} while (!flg && next_permutation(perms, perms + num));

	// If we found one
	if (flg) {
		*k = (num + m - 1) / m; // ceiling of num / m
		*phases = new char [*k * m];
		*basis = new boost::dynamic_bitset<> [*k * m];
		for (it = lst.begin(), i = 0; it != lst.end(); it++, i++) {
			(*phases)[perms[i]] = it->first;
			(*basis)[perms[i]]  = it->second;
		}

		// If the lst was not divisible by m
		if (*k != num / m) {
			for (int i = num; i < *k * m; i++) { 
				(*phases)[i] = 0;
			 	(*basis)[i]  = boost::dynamic_bitset<>(n);
			}
			make_full_rank(m, n, *basis + (*k - 1)*m);
		}

		delete [] tmp;
		return true;
	} else {
		delete [] tmp;
		return false;
	}
}

Circuit CNOTs_to_circuit(int m, int n, list<pair<int, int> > lst) {
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

Circuit compute_CNOT_network(int m, int n, const char * phases, const boost::dynamic_bitset<> * comp) {
	int i, j;
	// Compute steps to reduce to identity
	list<pair<int, int> > last  = gaussian(m, n, comp);

	// Do the inverse
	last.reverse();
	// Compute cnot network
	Circuit cnots = CNOTs_to_circuit(m, n, last);

	// Compute maximum number of T-stages
	int k = 0;
	for (i = 0; i < m; i++) k = max(k, (int)min((int)phases[i], 8 - (int)phases[i]));
	Circuit tees = Circuit(k);

	// Set the T gates
	for (i = 0; i < m; i++) {
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

Circuit compute_CNOT_network(int m, int n, const char * phases, 
		const boost::dynamic_bitset<> * comp, const boost::dynamic_bitset<> * inp) {
	int i, j;
	list<pair<int, int> > first = gaussian(m, n, inp);
	// Compute steps to reduce to identity
	list<pair<int, int> > last  = gaussian(m, n, comp);

	// Do the inverse
	last.reverse();
	first.splice(first.end(), last);
	// Compute cnot network
	Circuit cnots = CNOTs_to_circuit(m, n, first);

	// Compute maximum number of T-stages
	int k = 0;
	for (i = 0; i < m; i++) k = max(k, (int)min((int)phases[i], 8 - (int)phases[i]));
	Circuit tees = Circuit(k);

	// Set the T gates
	for (i = 0; i < m; i++) {
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

Circuit parallelize(int m, int n, int k, const char * phases, const boost::dynamic_bitset<> * bits) {
	Circuit acc(0), net, tmp;

	for (int i = 0; i < k; i++) {
		net = i == 0 ? compute_CNOT_network(m, n, phases, bits) 
								 : compute_CNOT_network(m, n, phases + i*m, bits + i*m, bits + (i-1)*m);
		tmp = acc;
		acc = tmp.append(net);
		delete_circuit(net);
		delete_circuit(tmp);
	}

	return acc;
}

/*
// Parallelize circuit INPUT on N qubits into a circuit with NUM_ANCILLA ancillas
Circuit T_parallelize(int n, const Circuit & input, int num_ancilla) {
	boost::dynamic_bitset<> * 

	// Parse the circuit into phase exponent form
	cout << "Parsing circuit... " << flush;
	cnot_t_circuit circ = parse_circuit(n, input);
	circ.print();


	cout << "Partitioning Matroid... " << flush;
	partitioning part = partition_matroid<
*/

int main() {
	num_qubits = 3;
	int n = 3;
	int m = 7;
	int k;
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

	cnot_t_circuit circ = parse_circuit(n, x);
	circ.print();

	/*
	char * phases;
	boost::dynamic_bitset<> * basis;

	cout << "Partitioning matroid\n" << flush;
	///////////////////////////////////////////////
	int num = ret.first.size();
	int h;
	boost::dynamic_bitset<> * teemp = new boost::dynamic_bitset<> [num];
	for (h = 0, it = ret.first.begin(); it != ret.first.end(); it++, h++) {
		teemp[h] = it->second;
	}
	partitioning liist = partition_matroid< boost::dynamic_bitset<> >(num, teemp, test_oracle);
	partitioning::iterator si;
	set<int>::iterator fi;
	for (si = liist.begin(); si != liist.end(); si++) {
		cout << "{";
		for (fi = si->begin(); fi != si->end(); fi++) {
			cout << *fi;
		}
		cout << "}";
	}
	cout << "\n";
	/////////////////////////////////////////////
	divide_into_bases(m, n, ret.first, &k, &phases, &basis);

	num_qubits = 7;
	cout << "Parallelizing\n" << flush;
	Circuit tmp = parallelize(m, n, k, phases, basis);

	// Add the final phase
	boost::dynamic_bitset<> * outperm = new boost::dynamic_bitset<> [m];
	for (int i = 0; i < m; i++) {
		if (i < n) {
			outperm[i] = ret.second[i];
		} else {
			outperm[i] = boost::dynamic_bitset<>(n);
		}
	}
	list<pair<int, int> > first = gaussian(m, n, basis + (k - 1)*m);
	// Compute steps to reduce to identity
	list<pair<int, int> > last  = gaussian(m, n, outperm);

	// Do the inverse
	last.reverse();
	first.splice(first.end(), last);
	// Compute cnot network
	Circuit cnots = CNOTs_to_circuit(m, n, first);
	Circuit final = tmp.append(cnots);

	final.print();

	*/
	return 1;
}
