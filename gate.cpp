#include "gate.h"

#define PHASE
/* Number of basis gates */
#define basis_size 9

/* Clifford group on 1 qubit */
#define I    0
#define H    1
#define X    2
#define Y    3
#define Z    4
#define S    5
#define Sd   6
#define T    7
#define Td   8

/* Basis state projections on one qubit */
#define PROJ(x, y)     (basis_size + 2*x + y)
#define IS_PROJ(x)  (basis_size <= x && x <= 0xfe)
#define GET_PROJ(x) (x - basis_size)

/* Controls. Highest order bit defines a control, the rest of the
   byte specifies the target */
#define C(x)          (x | 0x80)
#define IS_C(x)       (x & 0x80)
#define GET_TARGET(x) (x & 0x7F)

int num_qubits = 0;
int dim = 0;
int num_swaps = 0;

const string gate_names[] = {
  " I  ",
  " H  ",
  " X  ",  
  " Y  ",
  " Z  ",  
  " S  ",
  " S* ",
  " T  ",
  " T* ",
};

const char adjoint[] = {
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


Unitary * basis;
Unitary * swaps;

int fac(int n) {
  int ret = 1, i;
  for(i = n; i>1; i--) ret *= i;
  return ret;
}

int to_lexi(char * perm) {
  int i, tmp, ret = 0;
  for (i = 0 ; i < num_qubits-1; i++) {
    if (i == 0 || perm[i-1] <= perm[i]) {
      tmp = (perm[i] - i)*(num_qubits - 1 - i);
    } else {
      tmp = perm[i]*(num_qubits - 1 - i);
    }
    if (tmp < 0) tmp = 0;
    ret += tmp;
  }
  return ret;
}

char * from_lexi(int n) {
  char * ret = new char[num_qubits];
  int i, j, tmp, acc = n;
  for (i = 0; i < num_qubits; i++) {
    ret[i] = i;
  }
  for (i = 0; i < num_qubits-1; i++) {
    if (i == 0) {
      tmp = acc / (num_qubits - 1 - i);
    } else {
      tmp = acc / (num_qubits - 1 - i);
      if (tmp >= ret[i-1]) tmp += 1;
    }
    for (j = 0; j < num_qubits; j++) {
      if (ret[j] == tmp) {
        ret[j] = ret[i];
        ret[i] = tmp;
        break;
      }
    }
    acc = acc % (num_qubits - 1 - i);
  }
  return ret;
}


/* -------------- Gates */
Gate & Gate::operator=(const Gate & G) {
  for (int i = 0; i < num_qubits; i++) {
    gates[i] = G.gates[i];
  }
  return *this;
}

char & Gate::operator[](int i) const { 
  assert(0 <= i && i < num_qubits);
  return gates[i];
}

Gate::Gate() { gates = new char[num_qubits]; }
Gate::Gate(const Gate & G) { *this = G; }
Gate::~Gate() { delete[] gates; }

void Gate::adj(Gate & G) const {
  int i;

  for (i = 0; i < num_qubits; i++) {
    if (IS_C(gates[i])) {
      G[i] = gates[i];
    } else {
      G[i] = adjoint[gates[i]];
    }
  }
}

/* Return the tensor product of the 1-qubit matrices */
void Gate::tensor(Unitary & U) const {
  int i, j, k, x, y;
  LaComplex tmp;

  for (i = 0; i < num_qubits; i++) {
    assert(gates[i] < basis_size + dim);
  }

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      tmp = LaComplex(1, 0);
      x = i; y = j;
      for (k = num_qubits-1; k >= 0; k--) {
        tmp *= (LaComplex)(basis[gates[k]](x % 2, y % 2));
        x /= 2; y /= 2;
      }
      U(i, j) = tmp;
    }
  }
}

void Gate::print() const {
  int i;
  char tmp;

  for (i = 0; i < num_qubits; i++) {
    tmp = gates[i];

    if (IS_C(tmp)) {
      cout << "C(" << (int)GET_TARGET(tmp) + 1 << ")";
    } else {
      assert(tmp <= basis_size);
      cout << gate_names[tmp];
    }

    cout << "\n";
  }
}

/* Compute a unitary for the given gate */
void Gate::to_Unitary(Unitary & U) const {
  int i;
  for (i = 0; i < num_qubits; i++) {
    if (IS_C(gates[i])) {
      Gate A, B;
      Unitary V(dim, dim);
      char tmp = GET_TARGET(gates[i]);

      for (int j = 0; j < num_qubits; j++) {
        if (j == i) {
          A[j] = PROJ(i, 0);
          B[j] = PROJ(i, 1);
        } else if (j == tmp) {
          A[j] = I;
          B[j] = gates[j];
        } else {
          A[j] = gates[j];
          B[j] = gates[j];
        }
      }

      A.to_Unitary(V);
      B.to_Unitary(U);
      Blas_Mat_Mat_Mult(Unitary::eye(dim), V, U, false, false, 1, 1);
      return;
    }
  }
  (*this).tensor(U);
}

void Gate::permute(Gate & G, char * perm) const {
  for(int i = 0; i < num_qubits; i++) {
    G[i] = gates[perm[i]];
  }
}

/* --------------- Circuits */
void Circuit::print() const {
  const Circuit * tmp;
  char g;

  for (int i = 0; i < num_qubits; i++) {
    tmp = this;
    while (tmp) {
      g = tmp->G[i];
      if (IS_C(g)) {
        cout << "C(" << (int)GET_TARGET(g) + 1 << ")";
      } else {
        assert(g <= basis_size);
        cout << gate_names[g];
      }
      tmp = tmp->next;
    }
    cout << "\n";
  }
}

Circuit * Circuit::adj(Circuit * last) const {
  Circuit * tmp = new Circuit;
  tmp->next = last;
  G.adj(tmp->G);

  if (next) return next->adj(tmp);
  else      return tmp;
}

void Circuit::to_Unitary(Unitary & U) const {
  if (next != NULL) {
	  Unitary A(dim, dim);
    Unitary B(dim, dim);
    G.to_Unitary(A);
    next->to_Unitary(B);
    Blas_Mat_Mat_Mult(A, B, U, false, false, 1, 0);
  } else {
	  G.to_Unitary(U);
  }
}

/*---------------------------------*/

void swap_qubits(char * perm, Unitary & swap) {
  int i, j, acc, index;
  char * col = new char[num_qubits];
  for (i = 0; i < dim; i++) {
    index = 0;
    acc = i;
    for (j = 0; j < num_qubits; j++) {
      col[j] = acc / pow(2, num_qubits - 1 - j);
      acc %= (int)pow(2, num_qubits - 1 - j);
    }

    for (j = 0; j < num_qubits; j++) {
      index += col[perm[j]] * (pow(2, num_qubits - 1 - j));
    }

    swap(index, i) = LaComplex(1);
  }
}

void init(int n) {
  num_qubits = n;
  num_swaps = fac(n);
  dim = pow(2, n);
  basis = new Unitary[basis_size + dim];
  swaps = new Unitary[num_swaps];

  int i;
  double rt = 1.0/sqrt(2.0);

	for (i = 1; i < basis_size + dim; i++) {
	  basis[i] = Unitary::zeros(2);
	}

  basis[I] = Unitary::eye(2);

	basis[H](0, 0) = LaComplex(rt, 0);
	basis[H](1, 0) = LaComplex(rt, 0);
	basis[H](0, 1) = LaComplex(rt, 0);
	basis[H](1, 1) = LaComplex(-rt, 0);

	basis[X](0, 1) = LaComplex(1, 0);
	basis[X](1, 0) = LaComplex(1, 0);

	basis[Y](0, 1) = LaComplex(0, -1);
	basis[Y](1, 0) = LaComplex(0, 1);

	basis[Z](0, 0) = LaComplex(1, 0);
	basis[Z](1, 1) = LaComplex(-1, 0);

	basis[S](0, 0) = LaComplex(1, 0);
	basis[S](1, 1) = LaComplex(0, 1);

  Blas_Mat_Mat_Mult(basis[I],   basis[S],   basis[Sd], false, true, 1, 0);

	basis[T](0, 0) = LaComplex(1, 0);
	basis[T](1, 1) = LaComplex(rt, rt);

  Blas_Mat_Mat_Mult(basis[I],   basis[T],   basis[Td], false, true, 1, 0);

  for(i = 0; i < dim; i++) {
    basis[PROJ((i / 2), (i % 2))](i%2, i%2) = LaComplex(1, 0);
  }

  for(i = 0; i < num_swaps; i++) {
    swaps[i] = Unitary::zeros(dim);
    swap_qubits(from_lexi(i), swaps[i]);
    /*
    t1 = from_lehmer(i);
    t2 = to_lehmer(t1);
    cout << t2 << " " << i << ": ";
    for (int j = 0; j < num_qubits; j++) {
      cout << (int)t1[j] << ",";
    }
    cout << "\n";
    */
  }
}

double spec_norm(const Unitary & U) {
  Unitary G(U.size(0), U.size(1));
  Unitary F(U.size(0), U.size(1));
  LaVectorComplex e(G.size(1));
  double max = 0, temp;

  Blas_Mat_Mat_Mult(U, U, G, true, false, 1, 0);
  LaEigSolve(G, e, F);
 
  for (int i = 0; i < e.size(); i++) {
    temp = abs((complex<double>)((LaComplex)e(i)));
    temp = sqrt(temp);
    if (temp >= max) max = temp;
  }

  return max;
}

double dist(const Unitary & U, const Unitary & V) {
  #ifdef PHASE
    return spec_norm(U - V);
  #else
    double x, y, z;
    Unitary A(U.size(0), U.size(0));
    Unitary B(V.size(0), V.size(0));

    Blas_Mat_Mat_Mult(basis[X], U, A, false, true, 1, 0);
    Blas_Mat_Mat_Mult(U, A, A, false, false, 1, 0);
    Blas_Mat_Mat_Mult(basis[X], V, B, false, true, 1, 0);
    Blas_Mat_Mat_Mult(V, B, B, false, false, 1, 0);
    x = spec_norm(A - B);

    Blas_Mat_Mat_Mult(basis[Y], U, A, false, true, 1, 0);
    Blas_Mat_Mat_Mult(U, A, A, false, false, 1, 0);
    Blas_Mat_Mat_Mult(basis[Y], V, B, false, true, 1, 0);
    Blas_Mat_Mat_Mult(V, B, B, false, false, 1, 0);
    y = spec_norm(A - B);

    Blas_Mat_Mat_Mult(basis[Z], U, A, false, true, 1, 0);
    Blas_Mat_Mat_Mult(U, A, A, false, false, 1, 0);
    Blas_Mat_Mat_Mult(basis[Z], V, B, false, true, 1, 0);
    Blas_Mat_Mat_Mult(V, B, B, false, false, 1, 0);
    z = spec_norm(A - B);

    return sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
  #endif
}

int Hash_Unitary(const Unitary &U) {
  return (int)(10000000*dist(U, LaGenMatComplex::eye(dim)));
}

Canon canonicalize(const Unitary & U) {
  int i, d;
  int min = numeric_limits<int>::max();
  Unitary V(dim, dim), Vadj(dim, dim);

  Canon acc;

  for(i = 0; i < num_swaps; i++) {
    Blas_Mat_Mat_Mult(swaps[i], U, V, false, false, 1, 0);
    Blas_Mat_Mat_Mult(swaps[i], U, Vadj, false, true, 1, 0);

    d = Hash_Unitary(V);
    if (d < min) {
      min = d;
      acc.clear();
      acc.push_front(make_pair(0, i));
    } else if (d == min) {
      acc.push_front(make_pair(0, i));
    } 

    d = Hash_Unitary(Vadj);
    if (d < min) {
      min = d;
      acc.clear();
      acc.push_front(make_pair(1, i));
    } else if (d == min) {
      acc.push_front(make_pair(1, i));
    }
  }

  return acc;
} 

void test() {
  init(2);

  Circuit * A = new Circuit;
  Circuit * B = new Circuit;
  Circuit * C = new Circuit;
  Circuit * D = new Circuit;
  Circuit * E = new Circuit;
  Circuit * F = new Circuit;
  Circuit * G = new Circuit;
  A->G[0] = I;
  A->G[1] = S;
  B->G[0] = I;
  B->G[1] = H;
  C->G[0] = I;
  C->G[1] = T;
  D->G[0] = C(1);
  D->G[1] = X;
  E->G[0] = I;
  E->G[1] = Td;
  F->G[0] = I;
  F->G[1] = H;
  G->G[0] = I;
  G->G[1] = Sd;

  A->next = B;
  B->next = C;
  C->next = D;
  D->next = E;
  E->next = F;
  F->next = G;
  G->next = NULL;

  Unitary U(dim, dim);
  A->to_Unitary(U);
  A->print();
  cout << U << "\n";
  list< pair<char, int> > lst = canonicalize(U);
  list< pair<char, int> >::iterator iter = lst.begin();
  for(iter; iter != lst.end(); iter++) {
  pair<char, int> p = *iter;
  cout << (int)p.first << " " << (int)p.second << "\n";
  }

}
