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

void Rmatrix::resize(int mp, int np) {
  int i;

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
      mat[m][n] = M.mat[m][n];
    }
  }

  return *this;
}

Rmatrix & Rmatrix::operator+= (const Rmatrix & M) {
  int i, j;
  assert (m == M.m && n == M.n);

	for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      mat[m][n] += M.mat[m][n];
    }
  }

  return *this;
}

Rmatrix & Rmatrix::operator-= (const Rmatrix & M) {
  int i, j;
  assert (m == M.m && n == M.n);

	for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      mat[m][n] -= M.mat[m][n];
    }
  }

  return *this;
}

Rmatrix & Rmatrix::operator*= (const Elt & R) {
  int i, j;

	for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      mat[m][n] *= R;
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
    mat[i] = new Elt[n];
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

Elt & Rmatrix::operator() (int i, int j) {
  return mat[i][j];
}

Unitary Rmatrix::to_Unitary() const {
  int i, j;
  Unitary ret(m, n);

  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      ret(i, j) = LaComplex(mat[i][j].to_complex());
    }
  }

  return ret;
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
