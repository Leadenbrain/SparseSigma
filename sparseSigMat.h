// This is our actual sparse matrix.
// It will be an array of unordered maps
// Each member of the array will be a sparse vector representing each equation
// in the DAE.
#ifndef INCLUDED_SPARSESIGMAT_H
#define INCLUDED_SPARSESIGMAT_H

#include "sparseSigma.h"
/*
Class Constructor for the Sigma Matrix

Publically, we want:
    -Constructor: SigmaMatrix(int size, SigmaMatrixFcn dae_fun)
    -Destructor: ~SigmaMatrix()
    -Function to get sigma matrix: compSparseSigmaMatrix(SigmaMatrixFCN dae_fun
void* p) -> implemented in compSpSigma.cc
    -Function to print out sigma matrix printSparseSigmaMatrix() const;
        -> implemented in compSpSigma.cc

We will protect:
    -The DAE
    -Sigma Matrix
    -Size (n_)
*/
class SigmaMatrix {
 public:
  SigmaMatrix(int size, SigmaMatrixFcn dae_fun);

  void compSparseSigmaMatrix(SigmaMatrixFcn dae_fun, void* p);
  void printSparseSigmaMatrix() const;
  virtual ~SigmaMatrix();

 protected:
  SigmaMatrixFcn dae_fun_;
  SparseVecInt* sigma_matrix_;
  int n_;
};
#endif