#include "matrix.h"
#include <assert.h>

Rmatrix::Rmatrix() { m = n = 0; mat = NULL; }
Rmatrix::Rmatrix(int a, int b) { 
  m = a;
  n = b;
  mat = new Elt*[m];
  for (int i = 0; i < m; i++) {
    mat[i] = new Elt[n];
  }
}
Rmatrix::Rmatrix(const Rmatrix & M) { *this = M; }
Rmatrix::~Rmatrix() {
  for (int i = 0; i < m; i++) {
    delete [] mat[i];
  }
  delete [] mat;
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

Rmatrix &Rmatrix::operator*= (const Rmatrix & M) {
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
      for (k = 0; k < n; k++) {
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


Elt & Rmatrix::operator() (int i, int j) {
  return mat[i][j];
}

const bool Rmatrix::phase_eq(const Rmatrix & M) const {
  if (m != M.m || n != M.n) return false;
  Elt phase;

  int k, l, i, j, p;
  for (k = 0; k < 8; k++) {
    if (k/4 == 0) p = 1; else p = -1;
    if (k%4 == 0) phase = Elt(p, 0, 0, 0, 0);
    else if (k%4 == 1) phase = Elt(0, p, 0, 0, 0);
    else if (k%4 == 2) phase = Elt(0, 0, p, 0, 0);
    else phase = Elt(0, 0, 0, p, 0);

    for (l = 0; l < m*n; l++) {
      i = l / n;
      j = l % n;
      if (!(phase*mat[i][j] == M.mat[i][j])) break;
    }
    if (l == m*n) return true;
  }
  return false;
}

Unitary Rmatrix::to_Unitary() const {
  int i, j;
  Unitary ret(m, n);
  this->to_Unitary(ret);
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

void Rmatrix::adj(Rmatrix & M) const {
  int i, j;
  if (m != M.n || n != M.m) {
    M.resize(n, m);
  }
  for (i = 0; i < M.n; i++) {
    for (j = 0; j < M.m; j++) {
      M.mat[i][j] = mat[j][i].conj();
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
