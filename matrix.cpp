#include "matrix.h"
#include <assert.h>
#include <string.h>
#include <blas3pp.h>

/* Permutation related stuff ---------------------------------*/
const char * const * permutations;
const int * inversions;
const int * const * basis_permutations;

int fac(int n) {
  int ret = 1, i;
  for(i = n; i>1; i--) ret *= i;
  return ret;
}

const char * from_lexi(int i) {
  return permutations[(i < 0) ? inversions[abs(i)] : i];
}

int to_lexi(const char * perm) {
  int i, j;
  for (i = 0; i < num_perms; i++) {
    j = 0;
    while (perm[j] == permutations[i][j]) {
      j++;
      if (j >= num_qubits) return i;
    }
  }
  return -1;  
}

int permutation_helper(int mask, int ret, int index, char ** perms) {
  int i, tmp;
  for (i = 0; i < num_qubits; i++) {
    if ((1 << i) & mask) {
      if (index < (num_qubits - 1)) {
        tmp = permutation_helper(mask & ~(1 << i), ret, index + 1, perms);
        for (; ret < tmp; ret++) {
          perms[ret][index] = i;
        }
      } else {
        perms[ret][index] = i;
        ret += 1;
      }
    }
  }
  return ret;
}

/* Rmatrix stuff ----------------------------------------------*/
Rmatrix::Rmatrix() { m = n = 0; mat = NULL; }
Rmatrix::Rmatrix(int a, int b) { 
  m = a;
  n = b;
  mat = new Elt[m*n];
}
Rmatrix::Rmatrix(const Rmatrix & M) {
  int i, j;

  m = M.m; 
  n = M.n; 
  mat = new Elt[m*n];
  for (i = 0; i < m*n; i++) {
    mat[i] = M.mat[i];
  }
}
  
Rmatrix::~Rmatrix() {
  delete [] mat;
}

Rmatrix Rmatrix::rand(int m, int n) {
  int i, j;
  Rmatrix ret(m, n);
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      ret(i, j) = Elt::randelt();
      ret(i, j).reduce();
    }
  }
  return ret;
}

Rmatrix zero(int m, int n) {
  int i, j;
  Rmatrix ret(m, n);
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      ret(i, j) = Elt(0, 0, 0, 0, 0);
    }
  }
  return ret;
}

Rmatrix eye(int m, int n) {
  int i, j;
  Rmatrix ret(m, n);
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      if (i == j) ret(i, j) = Elt(1, 0, 0, 0, 0);
      else ret(i, j) = Elt(0, 0, 0, 0, 0);
    }
  }
  return ret;
}

Rmatrix col_permutation(int m, int n, int perm) {
  Rmatrix ret = zero(m, n);
  for (int i = 0; i < m; i++) {
    ret(i, basis_permutations[perm][i]) = Elt(1, 0, 0, 0, 0);
  }
  return ret;
}

int Rmatrix::rows() const {
  return m;
}

int Rmatrix::cols() const {
  return n;
}

void Rmatrix::resize(int mp, int np) {
  int i;

  delete [] mat;

  m = mp;
	n = np;
	mat = new Elt[m*n];
}

Rmatrix & Rmatrix::operator= (const Rmatrix & M) {
  int i, j;
	if (m != M.m || n != M.n) {
    this->resize(M.m, M.n);
  }

	for (i = 0; i < m*n; i++) {
    mat[i] = M.mat[i];
  }

  return *this;
}

Rmatrix & Rmatrix::operator+= (const Rmatrix & M) {
  int i, j;
  assert (m == M.m && n == M.n);

	for (i = 0; i < m*n; i++) {
    mat[i] += M.mat[i];
    mat[i].reduce();
  }

  return *this;
}

Rmatrix & Rmatrix::operator-= (const Rmatrix & M) {
  int i, j;
  assert (m == M.m && n == M.n);

	for (i = 0; i < m*n; i++) {
    mat[i] -= M.mat[i];
    mat[i].reduce();
  }

  return *this;
}

Rmatrix & Rmatrix::operator*= (const Elt & R) {
  int i, j;

	for (i = 0; i < m*n; i++) {
    mat[i] *= R;
    mat[i].reduce();
  }

  return *this;
}

/* unfinished Strassen algorithm
Rmatrix Rmatrix::FMM(const Rmatrix & M) {
  for (int k = m >> 1; k > 0; k = k >> 1) {
    for (int i = 0; i < k; i++) {
      for (int j = 0; j < k; j++) {
        M1 = mat[


  if (m == 2) {
    Elt M1 = (mat[0][0]+mat[1][1])*(M.mat[0][0]+M.mat[1][1]);
    Elt M2 = (mat[1][0]+mat[1][1])*(M.mat[0][0]);
    Elt M3 = (mat[0][0])*(M.mat[0][1]+M.mat[1][1]);
    Elt M4 = (mat[1][1])*(M.mat[1][0]-M.mat[0][0]);
    Elt M5 = (mat[0][0]+mat[0][1])*(M.mat[1][1]);
    Elt M6 = (mat[1][0]-mat[0][0])*(M.mat[0][0]+M.mat[0][1]);
    Elt M7 = (mat[0][1]-mat[1][1])*(M.mat[1][0]+M.mat[1][1]);
    Rmatrix ret(2, 2);
    ret.mat[0][0] = M1+M4-M5+M7;
    ret.mat[0][1] = M3+M5;
    ret.mat[1][0] = M2+

}
*/

Rmatrix & Rmatrix::operator*= (const Rmatrix & M) {
  int i, j, k;
  Elt sum;
  assert (n == M.m);

  if (n != M.n) {
    // Allocates new memory
    Elt * newmat;

    newmat = new Elt[m*M.n];

    for (j = 0; j < M.n; j++) {
      for (i = 0; i < m; i++) {
        sum = Elt(0, 0, 0, 0, 0);
        for (k = 0; k < M.m; k++) {
          sum += mat[i*n + k]*M.mat[k*M.n + j];
          sum.reduce();
        }
        newmat[i*M.n + j] = sum;
      }
    }

    if (n != M.n) n = M.n;
    delete [] mat;

    mat = newmat;
    return *this;
  } else {
    // Faster, does no reallocate
    Elt tmp[n];

    for (i = 0; i < m; i++) {
      for (j = 0; j < M.n; j++) {
        sum = Elt(0, 0, 0, 0, 0);
        for (k = 0; k < M.m; k++) {
          sum += mat[i*n + k]*M.mat[k*M.n + j];
          sum.reduce();
        }
        tmp[j] = sum;
      }
      for (j = 0; j < M.n; j++) {
        mat[i*n + j] = tmp[j];
      }
    }
    return *this;
  }
}

Rmatrix &Rmatrix::left_multiply(const Rmatrix & M) {
  int i, j, k;
  Elt sum;
  assert (n == M.m);
  if (n != M.n) {
    // Allocates new memory
    Elt * newmat;

    newmat = new Elt[M.m*n];

    for (j = 0; j < n; j++) {
      for (i = 0; i < M.m; i++) {
        sum = Elt(0, 0, 0, 0, 0);
        for (k = 0; k < m; k++) {
          sum += M.mat[i*M.n+k]*mat[k*n+j];
          sum.reduce();
        }
        newmat[i*n+j] = sum;
      }
    }

    if (m != M.m) m = M.m;
    delete [] mat;

    mat = newmat;
    return *this;
  } else {
    // Faster, does no reallocate
    Elt tmp[n];

    for (j = 0; j < n; j++) {
      for (i = 0; i < M.m; i++) {
        sum = Elt(0, 0, 0, 0, 0);
        for (k = 0; k < m; k++) {
          sum += M.mat[i*M.n+k]*mat[k*n+j];
          sum.reduce();
        }
        tmp[i] = sum;
      }
      for (i = 0; i < m; i++) {
        mat[i*n+j] = tmp[i];
      }
    }
    return *this;
  }
}

const Rmatrix Rmatrix::operator+ (const Rmatrix & M) const {
  Rmatrix ret = *this;
  ret += M;
  return ret;
}
const Rmatrix Rmatrix::operator- (const Rmatrix & M) const {
  Rmatrix ret = *this;
  ret -= M;
  return ret;
}
const Rmatrix Rmatrix::operator* (const Elt & R) const {
  Rmatrix ret = *this;
  ret *= R;
  return ret;
}
const Rmatrix Rmatrix::operator* (const Rmatrix & M) const {
  Rmatrix ret = *this;
  ret *= M;
  return ret;
}

const bool Rmatrix::operator== (const Rmatrix & M) const {
  if (m != M.m || n != M.n) return false;

  int i, j;
  for (i = 0; i < m*n; i++) {
    if (not (mat[i] == M.mat[i])) return false;
  }
  return true;
}

const bool Rmatrix::operator<  (const Rmatrix & M) const {
  if (m != M.m || n != M.n) {
    if (m < M.m) return true;
    if (m > M.m) return false;
    if (n < M.n) return true;
    if (n > M.n) return false;
  }

  int i, j;
  for (i = 0; i < m*n; i++) {
    if (mat[i] < M.mat[i]) return true;
    if (not (mat[i] == M.mat[i])) return false;
  }
  return false;
}

const Elt & Rmatrix::operator() (int i, int j) const {
	if (mat == NULL) {
		cout << "ERROR: unallocated matrix\n";
		exit(1);
	}
  return mat[i*n+j];
}

Elt & Rmatrix::operator() (int i, int j) {
	if (mat == NULL) {
		cout << "ERROR: unallocated matrix\n";
		exit(1);
	}
  return mat[i*n+j];
}

Elt Rmatrix::trace() const {
	Elt ret(0, 0, 0, 0, 0);

	if (m != n) {
		cout << "ERROR: trace of non-square matrix\n";
		exit(1);
	}
	for (int i = 0; i < m; i++) {
		ret += mat[i*n + i];
	}
	ret.reduce();
	return ret;
}

const bool Rmatrix::phase_eq(const Rmatrix & M) const {
  if (m != M.m || n != M.n) return false;
  Elt phase(0, 1, 0, 0, 0), ph(0, 1, 0, 0, 0);

  int k, l, i, j, p;
  bool flg;
  /*
  for (k = 0; k < 8; k++) {
    if (k/4 == 0) p = 1; else p = -1;
    if (k%4 == 0) phase = Elt(p, 0, 0, 0, 0);
    else if (k%4 == 1) phase = Elt(0, p, 0, 0, 0);
    else if (k%4 == 2) phase = Elt(0, 0, p, 0, 0);
    else phase = Elt(0, 0, 0, p, 0);

    flg = true;
    for (i = 0; i < m && flg; i++) {
      for (j = 0; j < n && flg; j++) {
        flg = flg && (phase*mat[i][j] == M.mat[i][j]);
      }
    }
    if (flg == true) return true;
  }
  return false;
  */
  k = 0;
  for (i = 0; i < m*n; i++) {
    if (mat[i].is_zero()) {
      if (M.mat[i].is_zero()) {
        break;
      } else {
        return false;
      }
    } else {
      while (!(ph*mat[i] == M.mat[i])) {
        k++;
        if (k >= 8) return false;
        ph *= phase;
      }
    }
  }
  return true;
}

const bool Rmatrix::equiv(const Rmatrix & M) const {
  int i;
  Rmatrix t(dim, dim);
  for (i = 0; i < 2*num_perms; i++) {
    if (i % 2 == 0) {
      M.permute(t, i / 2);
    } else {
      M.permute_adj(t, i / 2);
    }
    if (this->phase_eq(t)) return true;
  }
  return false;
}

Unitary Rmatrix::to_Unitary() const {
  Unitary ret(m, n);
  this->to_Unitary(ret);
  return ret;
}

LaGenMatDouble Rmatrix::to_Unitary_abs() const {
  LaGenMatDouble ret(m, n);
  this->to_Unitary_abs(ret);
  return ret;
}

Unitary Rmatrix::to_Unitary_canon() const {
  Unitary ret(m, n);
  this->to_Unitary_canon(ret);
  return ret;
}

void Rmatrix::to_Unitary(Unitary & U) const {
  int i, j;
  if(U.size(0) != m || U.size(1) != n) {
    U.resize(m, n);
    cout << "Resize\n";
  }
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      U(i, j) = LaComplex(mat[i*n+j].to_complex());
    }
  }
}

void Rmatrix::to_Unitary_abs(LaGenMatDouble & U) const {
  int i, j;
  if(U.size(0) != m || U.size(1) != n) {
    U.resize(m, n);
    cout << "Resize\n";
  }
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      U(i, j) = mat[i*n+j].abs();
    }
  }
}

void Rmatrix::to_Unitary_canon(Unitary & U) const {
  int i, j;
  bool flg = false;
  Elt phase, tmp;
  if(U.size(0) != m || U.size(1) != n) {
    U.resize(m, n);
    cout << "Resize\n";
  }
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      if (mat[i*n+j].is_zero()) {
        U(i, j) = LaComplex(0, 0);
      } else if (!flg) {
        flg = true;
        phase = mat[i*n+j].conj();
        tmp = phase*mat[i*n+j];
        tmp.reduce();
        U(i, j) = LaComplex(tmp.to_complex());
      } else {
        tmp = phase*mat[i*n+j];
        tmp.reduce();
        U(i, j) = LaComplex(tmp.to_complex());
      }
    }
  }
}

void Rmatrix::adj(Rmatrix & M) const {
  int i, j;
  if (M.m != n || M.n != m) {
    M.resize(n, m);
  }
  for (i = 0; i < M.m; i++) {
    for (j = 0; j < M.n; j++) {
      M.mat[i*M.n+j] = mat[j*n+i].conj();
    }
  }
}

void Rmatrix::permute(Rmatrix & M, int x) const {
  int i, j, ip, jp;
  if (M.m != m || M.n != n) {
    M.resize(m, n);
  }
  for (i = 0; i < m; i++) {
    ip = basis_permutations[x][i];
    for (j = 0; j < n; j++) {
      jp = basis_permutations[x][j];
      if (ip < m && jp < n) {
        M.mat[ip*M.n+jp] = mat[i*n+j];
      }
    }
  }
}

void Rmatrix::permute_adj(Rmatrix & M, int x) const {
  int i, j, ip, jp;
  if (M.m != n || M.n != m) {
    M.resize(n, m);
  }
  for (i = 0; i < m; i++) {
    ip = basis_permutations[x][i];
    for (j = 0; j < n; j++) {
      jp = basis_permutations[x][j];
      M.mat[ip*M.n+jp] = mat[j*n+i].conj();
    }
  }
}

void Rmatrix::print() const {
  int i, j;
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      mat[i*n+j].print();
      cout << " ";
    }
    cout << "\n";
  }
}

void Rmatrix::submatrix(int x, int y, int numrow, int numcol, Rmatrix & M) const {
  int i, j;
  if (M.m != numrow || M.n != numcol) {
    M.resize(numrow, numcol);
  }
  for (i = 0; i < numrow; i++) {
    for (j = 0; j < numcol; j++) {
      M.mat[i*M.n+j] = mat[(i+x)*n+j+y];
    }
  }
}

void Rmatrix::canon_phase() {
  Elt phase;
  bool flg = false;
  for (int i = 0; i < m*n; i++) {
    if (!mat[i].is_zero()) {
      if (!flg) {
        flg = true;
        phase = mat[i].conj();
      }
      mat[i] *= phase;
      mat[i].reduce();
    }
  }
}

inline bool parity(unsigned char b) {
  b ^= b>4;
  b ^= b>2;
  b ^= b>1;
  return b & 1;
}

inline unsigned char apply_fnc(unsigned char * fnc, unsigned char x, int bits) {
  int i;
  unsigned char ret = 0;
  for (i = 0; i < bits; i++) {
    ret = (ret << 1) | parity(fnc[i] & x);
  }
  return ret;
}


bool Rmatrix::is_nonlinear_reversible() const {
  int i, j;
  unsigned char fnc[m], tmp1, tmp2;
  /*
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      if (mat[i][j] == Elt(0, 0, 0, 0, 0)) {
        fnc[i] = fnc[i] << 1;
      } else if (mat[i][j] == Elt(1, 0, 0, 0, 0)) {
        fnc[i] = (fnc[i] << 1) | 1;
      } else {
        return false;
      }
    }
  }

  for (i = 0; i < (1 << n) - 1; i++) {
    for (j = i+1; j < (1 << n); j++) {
      cout << "HELLO " << i << " " << j << "\n";
      tmp1 = apply_fnc(fnc, i, m) ^ apply_fnc(fnc, j, m);
      tmp2 = apply_fnc(fnc, i ^ j, m);
      if (tmp1 != tmp2) return true;
    }
  }
  */
  return false;
}

void init_rmatrix() {
  int i, j;

  /* Allocate some memory... */
  char ** permutations_tmp = new char*[num_perms];
  int * inversions_tmp = new int[num_perms];
  int ** basis_permutations_tmp = new int*[num_perms];
  
  for (i = 0; i < num_perms; i++) {
    permutations_tmp[i] = new char[num_qubits];
    basis_permutations_tmp[i] = new int[1 << num_qubits];
  }

  /* Generate permutations */
  permutation_helper(~((int)0), 0, 0, permutations_tmp);
  permutations = permutations_tmp;

  /* Generate inversions */
  char tmp[num_qubits];
  for (i = 0; i < num_perms; i++) {
    for (j = 0; j < num_qubits; j++) {
      tmp[permutations[i][j]] = j;
    }
    inversions_tmp[i] = to_lexi(tmp);
  }
  inversions = inversions_tmp;

  /* Generate basis state permutations */
  int tmp2, k;
  for (i = 0; i < num_perms; i++) {
    for (j = 0; j < (1 << num_qubits); j++) {
      tmp2 = 0;
      /* To binary, permutated */
      for (k = 0; k < num_qubits; k++) {
        tmp[permutations[i][k]] = (bool)((j << k) & (1 << (num_qubits - 1)));
      }
      for (k = 0; k < num_qubits; k++) {
        tmp2 |= tmp[k] << (num_qubits - 1 - k);
      }
      basis_permutations_tmp[i][j] = tmp2;
    }
  }
  basis_permutations = basis_permutations_tmp;
}


void test_rmatrix() {
  Rmatrix A;
  Rmatrix B(4, 4);
  Rmatrix C = eye(2, 2);
  Rmatrix D(C);
  assert(C == D);
  assert(C.rows() == 2);
  C.resize(6, 2);
  assert(C.rows() == 6);
  C = eye(6, 2);
  Rmatrix *X = new Rmatrix(2, 2);
  (*X)(0, 0) = Elt(0, 1, 0, 0, 0);
  (*X)(0, 1) = Elt(0, 1, 0, 0, 0);
  (*X)(1, 0) = Elt(0, 1, 0, 0, 0);
  (*X)(1, 1) = Elt(0, 1, 0, 0, 0);
  assert(D*(*X) == (*X));
  A = C*(*X);
  assert(A.rows() == 6);
  assert(X->phase_eq((*X)*Elt(0, 0, 1, 0, 0)));
  X->adj(B);
  assert(X->phase_eq(B));
  X->canon_phase();
  B.canon_phase();
  assert(B == *X);
  /*
  A.submatrix(0, 0, 2, 2, B);
  assert(B == *X);
  C *= (*X);
  assert(A == C);
*/
  delete X;
}

//----------------------------------------------------------------

void adj_unitary(const Unitary & A, Unitary & B) {
  int i, j;
  if (B.rows() != A.cols() || B.cols() != A.rows()) {
    B.resize(A.cols(), A.rows());
  }
	Blas_Mat_Mat_Mult(Unitary::eye(A.rows()), A, B, false, true, 1, 0);
}

void permute_unitary(const Unitary & A, Unitary & B, int x) {
  int i, j, ip, jp;
  if (B.rows() != A.cols() || B.cols() != A.rows()) {
    B.resize(A.cols(), A.rows());
  }
  for (i = 0; i < A.rows(); i++) {
    ip = basis_permutations[x][i];
    for (j = 0; j < A.cols(); j++) {
      jp = basis_permutations[x][j];
      if (ip < A.rows() && jp < A.cols()) {
				B(ip, jp) = A(i, j);
      }
    }
  }
}
