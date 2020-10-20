// Quick and easy implementation of the matrices
// Computation contained in compSpSigma.cc
#include "sparseSigMat.h"
#include <stdexcept>

/*
Sparse Sigma Matrix Constructor

-Reads in a DAE and initializes an array of sparse matrices
-Initializes DAE as a private variable
*/
SigmaMatrix::SigmaMatrix(int size, SigmaMatrixFcn dae_fun) : n_(size) {
  sigma_matrix_ = new SparseVecInt[n_];

  dae_fun_ = dae_fun;
}

/*
Sparse Sigma Matrix Destructor

-Deletes the entries in the sparse matrix
-Deletes sparse matrix
*/
SigmaMatrix::~SigmaMatrix() {
  delete[] sigma_matrix_;
}