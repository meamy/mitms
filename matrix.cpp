#include "matrix.h"
#include <assert.h>
#include <string.h>

/* Permutation related stuff ---------------------------------*/
int num_elts;
int num_permutations;
char ** permutations;
int * inversions;
int ** basis_permutations;

int fac(int n) {
  int ret = 1, i;
  for(i = n; i>1; i--) ret *= i;
  return ret;
}

char * from_lexi_(int i) {
  return permutations[i];
}

int to_lexi_(char * perm) {
  int i, j;
  for (i = 0; i < num_permutations; i++) {
    j = 0;
    while (perm[j] == permutations[i][j]) {
      j++;
      if (j >= num_elts) return i;
    }
  }
  return -1;  
}

int permutation_helper(int mask, int ret, int index) {
  int i, tmp;
  for (i = 0; i < num_elts; i++) {
    if ((1 << i) & mask) {
      if (index < (num_elts - 1)) {
        tmp = permutation_helper(mask & ~(1 << i), ret, index + 1);
        for (; ret < tmp; ret++) {
          permutations[ret][index] = i;
        }
      } else {
        permutations[ret][index] = i;
        ret += 1;
      }
    }
  }
  return ret;
}

void init_permutations(int num) {
  int i, j;
  /* Allocate some memory... */
  num_elts = num;
  num_permutations = fac(num);
  permutations = new char*[num_permutations];
  inversions = new int[num_permutations];
  basis_permutations = new int*[num_permutations];
  for (i = 0; i < num_permutations; i++) {
    permutations[i] = new char[num_elts];
    basis_permutations[i] = new int[1 << num_elts];
  }

  /* Generate permutations */
  permutation_helper(~((int)0), 0, 0);
/*
  cout << "Permutations\n";
  for (i = 0; i < num_permutations; i++) {
    cout << i << ": ";
    for (j = 0; j < num_elts; j++) {
      cout << (int)permutations[i][j] << " ";
    }
    cout << "\n";
  }
  cout << "\n";
*/

  /* Generate inversions */
  char tmp[num];
//cout << "Inversion\n";
  for (i = 0; i < num_permutations; i++) {
    for (j = 0; j < num_elts; j++) {
      tmp[permutations[i][j]] = j;
    }
    inversions[i] = to_lexi_(tmp);
//  cout << i << ": " << inversions[i] << "\n";
  }
//cout << "\n";

  /* Generate basis state permutations */
  int tmp2, k;
//cout << "Basis states\n";
  for (i = 0; i < num_permutations; i++) {
//  cout << i << ": ";
    for (j = 0; j < (1 << num_elts); j++) {
      tmp2 = 0;
//    cout << j << "->";
      /* To binary, permutated */
      for (k = 0; k < num_elts; k++) {
        tmp[permutations[i][k]] = (bool)((j << k) & (1 << (num_elts - 1)));
      }
      for (k = 0; k < num_elts; k++) {
//      cout << (int)tmp[k];
        tmp2 |= tmp[k] << (num_elts - 1 - k);
      }
//    cout << "->" << tmp2 << " ";
      basis_permutations[i][j] = tmp2;
    }
//  cout << "\n";
  }
}

/* Rmatrix stuff ----------------------------------------------*/
Rmatrix::Rmatrix() { m = n = 0; mat = NULL; }
Rmatrix::Rmatrix(int a, int b) { 
  m = a;
  n = b;
  mat = new Elt*[m];
  for (int i = 0; i < m; i++) {
    mat[i] = new Elt[n];
  }
}
Rmatrix::Rmatrix(const Rmatrix & M) {
  int i, j;

  m = M.m; 
  n = M.n; 
  mat = new Elt*[m];
  for (i = 0; i < m; i++) {
    mat[i] = new Elt[n];
    for (j = 0; j < n; j++) {
      mat[i][j] = M.mat[i][j];
    }
  }
}
  
Rmatrix::~Rmatrix() {
  for (int i = 0; i < m; i++) {
    delete [] mat[i];
  }
  delete [] mat;
}

Rmatrix Rmatrix::rand(int m, int n) {
  int i, j;
  Rmatrix ret(m, n);
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      ret(i, j) = Elt::rand();
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

int Rmatrix::rows() const {
  return m;
}

int Rmatrix::cols() const {
  return n;
}

void Rmatrix::resize(int mp, int np) {
  int i;

  for (int i = 0; i < m; i++) {
    delete [] mat[i];
  }
  delete [] mat;

  m = mp;
	n = np;
	mat = new Elt*[m];
  for (i = 0; i < m; i++) {
    mat[i] = new Elt[n];
  }
}

Rmatrix & Rmatrix::operator= (const Rmatrix & M) {
  int i, j;
	if (m != M.m || n != M.n) {
    this->resize(M.m, M.n);
  }

	for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      mat[i][j] = M.mat[i][j];
    }
  }

  return *this;
}

Rmatrix & Rmatrix::operator+= (const Rmatrix & M) {
  int i, j;
  assert (m == M.m && n == M.n);

	for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      mat[i][j] += M.mat[i][j];
    }
  }

  return *this;
}

Rmatrix & Rmatrix::operator-= (const Rmatrix & M) {
  int i, j;
  assert (m == M.m && n == M.n);

	for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      mat[i][j] -= M.mat[i][j];
    }
  }

  return *this;
}

Rmatrix & Rmatrix::operator*= (const Elt & R) {
  int i, j;

	for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      mat[i][j] *= R;
    }
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

Rmatrix &Rmatrix::operator*= (const Rmatrix & M) {
  /*
  int i, j, k;
  Elt ** newmat, sum;
  assert (n == M.m);

  if (n != M.n) n = M.n;
	newmat = new Elt*[m];
  for (i = 0; i < m; i++) {
    newmat[i] = new Elt[n];
  }

  for (j = 0; j < M.n; j++) {
    for (i = 0; i < m; i++) {
      sum = Elt(0, 0, 0, 0, 0);
      for (k = 0; k < M.m; k++) {
        sum += mat[i][k]*M.mat[k][j];
      }
      newmat[i][j] = sum;
    }
  }

  for (int i = 0; i < m; i++) {
    delete [] mat[i];
  }
  delete [] mat;

  mat = newmat;
  return *this;
  */
  assert (n == M.m);
  assert (n == M.n);
  int i, j, k;
  Elt sum, tmp[n];

  for (i = 0; i < m; i++) {
    for (j = 0; j < M.n; j++) {
      sum = Elt(0, 0, 0, 0, 0);
      for (k = 0; k < M.m; k++) {
        sum += mat[i][k]*M.mat[k][j];
      }
      tmp[j] = sum;
    }
    for (j = 0; j < M.n; j++) {
      mat[i][j] = tmp[j];
    }
  }
  return *this;
}

Rmatrix &Rmatrix::left_multiply(const Rmatrix & M) {
  /*
  int i, j, k;
  Elt ** newmat, sum;
  assert (m == M.n);

  if (m != M.m) m = M.m;
	newmat = new Elt*[m];
  for (i = 0; i < m; i++) {
    newmat[i] = new Elt[n];
  }

  for (j = 0; j < M.n; j++) {
    for (i = 0; i < m; i++) {
      sum = Elt(0, 0, 0, 0, 0);
      for (k = 0; k < M.m; k++) {
        sum += M.mat[i][k]*mat[k][j];
      }
      newmat[i][j] = sum;
    }
  }

  for (int i = 0; i < m; i++) {
    delete [] mat[i];
  }
  delete [] mat;

  mat = newmat;
  return *this;
  */
  assert (n == M.m);
  assert (n == M.n);
  int i, j, k;
  Elt sum, tmp[n];

  for (j = 0; j < n; j++) {
    for (i = 0; i < M.m; i++) {
      sum = Elt(0, 0, 0, 0, 0);
      for (k = 0; k < m; k++) {
        sum += M.mat[i][k]*mat[k][j];
      }
      tmp[i] = sum;
    }
    for (i = 0; i < m; i++) {
      mat[i][j] = tmp[i];
    }
  }
  return *this;
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
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      if (not (mat[i][j] == M.mat[i][j])) return false;
    }
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
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      if (mat[i][j] < M.mat[i][j]) return true;
      if (not (mat[i][j] == M.mat[i][j])) return false;
    }
  }
  return false;
}

const Elt & Rmatrix::operator() (int i, int j) const {
  return mat[i][j];
}

Elt & Rmatrix::operator() (int i, int j) {
  return mat[i][j];
}

const bool Rmatrix::phase_eq(const Rmatrix & M) const {
  if (m != M.m || n != M.n) return false;
  Elt phase(0, 1, 0, 0, 0);

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
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      if (mat[i][j].is_zero()) {
        if (M.mat[i][j].is_zero()) {
          break;
        } else {
          return false;
        }
      } else {
        while (!(phase*mat[i][j] == M.mat[i][j])) {
          k++;
          if (k >= 8) return false;
          phase = phase * phase;
        }
      }
    }
  }
  return true;
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
      U(i, j) = LaComplex(mat[i][j].to_complex());
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
      U(i, j) = mat[i][j].abs();
    }
  }
}

void Rmatrix::to_Unitary_canon(Unitary & U) const {
  int i, j;
  bool flg = false;
  Elt phase;
  if(U.size(0) != m || U.size(1) != n) {
    U.resize(m, n);
    cout << "Resize\n";
  }
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      if (mat[i][j].is_zero()) {
        U(i, j) = LaComplex(0, 0);
      } else if (!flg) {
        flg = true;
        phase = mat[i][j].conj();
        U(i, j) = LaComplex((phase*mat[i][j]).to_complex());
      } else {
        U(i, j) = LaComplex((phase*mat[i][j]).to_complex());
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
      M.mat[i][j] = mat[j][i].conj();
    }
  }
}

void Rmatrix::permute(Rmatrix & M, int x) const {
  int i, j, ip, jp;
  if (M.m != m || M.n != n) {
    M.resize(n, m);
  }
  for (i = 0; i < m; i++) {
    ip = basis_permutations[x][i];
    for (j = 0; j < n; j++) {
      jp = basis_permutations[x][j];
      M.mat[ip][jp] = mat[i][j];
    }
  }
}

void Rmatrix::permute_adj(Rmatrix & M, int x) const {
  int i, j, ip, jp;
  if (M.m != m || M.n != n) {
    M.resize(n, m);
  }
  for (i = 0; i < m; i++) {
    ip = basis_permutations[x][i];
    for (j = 0; j < n; j++) {
      jp = basis_permutations[x][j];
      M.mat[ip][jp] = mat[j][i].conj();
    }
  }
}

void Rmatrix::print() const {
  int i, j;
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      mat[i][j].print();
      cout << " ";
    }
    cout << "\n";
  }
}

void Rmatrix::submatrix(int m, int n, int numrow, int numcol, Rmatrix & M) const {
  int i, j;
  if (M.m != numrow || M.n != numcol) {
    M.resize(numrow, numcol);
  }
  for (i = 0; i < numrow; i++) {
    for (j = 0; j < numcol; j++) {
      M.mat[i][j] = mat[i+m][j+n];
    }
  }
}

void Rmatrix::canon_phase() {
  Elt phase;
  bool flg = false;
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      if (!mat[i][j].is_zero()) {
        if (!flg) {
          flg = true;
          phase = mat[i][j].conj();
        }
        mat[i][j] *= phase;
      }
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

  for (i = 0; i < pow(2, n) - 1; i++) {
    for (j = i+1; j < pow(2, n); j++) {
      cout << "HELLO " << i << " " << j << "\n";
      tmp1 = apply_fnc(fnc, i, m) ^ apply_fnc(fnc, j, m);
      tmp2 = apply_fnc(fnc, i ^ j, m);
      if (tmp1 != tmp2) return true;
    }
  }
  return false;
}

void matrix_test() {
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
