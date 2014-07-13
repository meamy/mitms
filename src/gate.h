
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

#ifndef GATE
#define GATE

#include "matrix.h"
#include <cstring>

/* -------------- Gates */
typedef char* gate;

void next_gate(gate G);
void next_clifford(gate G);

bool gate_eq(const gate A, const gate B);
inline bool gate_neq(const gate A, const gate B) { return !(A == B); }
bool is_eye(const gate G);
bool nontrivial_id(const gate A, const gate B);

void gate_transform(const gate A, gate B, const char * perm, bool adj);
void gate_transform(const gate A, gate B, int i, bool adj);
void gate_to_Rmatrix(const gate G, Rmatrix & U, bool adj);
void gate_to_Unitary(const gate G, Unitary & U, bool adj);

void print_gate(const gate G);
void output_gate(const gate G, ofstream & out);
void input_gate(gate G, ifstream & in);

unsigned int gate_hasher(const gate R);

inline void copy_gate(const gate A, gate B)
  { memcpy(B, A, num_qubits); }
inline void gate_adj(const gate A, gate B)                         
  { gate_transform(A, B, 0, true); }
inline void gate_permute(const gate A, gate B, const char *  perm) 
  { gate_transform(A, B, perm, false); }
inline void gate_permute(const gate A, gate B, int i)              
  { gate_transform(A, B, i, false); }
inline void gate_permute_adj(const gate A, gate B, const char * perm) 
  { gate_transform(A, B, perm, true); }
inline void gate_permute_adj(const gate A, gate B, int i)
  { gate_transform(A, B, i, true); }
inline void gate_to_Rmatrix(const gate G, Rmatrix & U)
  { gate_to_Rmatrix(G, U, false); }
inline void gate_to_Unitary(const gate G, Unitary & U)
  { gate_to_Unitary(G, U, false); }

void init_gate();
void test_gate();

#endif
