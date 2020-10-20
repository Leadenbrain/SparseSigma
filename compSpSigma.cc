// A separate file for the implementation of the computation
#include "sparseSigMat.h"

/*
This function will actually compute our sigma matrix
  -Reads in the DAE and the parameters
  -Initializes and allocates vectors
  -Tries running the DAE to ensure no errors
  -Compute the Sigma Matrix as entries of the DAE functions
  -Delete f and y
*/
void SigmaMatrix::compSparseSigmaMatrix(SigmaMatrixFcn dae_fun,
                                        void* dae_params) {
  SigmaVector *y, *f;

  // Allocate y and f pointer, call default constructor.
  y = new SigmaVector[n_];
  f = new SigmaVector[n_];
  // initialize the SigmaVector vector y
  for (int i = 0; i < n_; i++)
    y[i].initialize_vector(i);

  int nnz = 0;
  // compute the SigmaVector matrix; check for errors
  try {
    dae_fun(y, f, dae_params);
  } catch (const std::exception& exc) {
    std::cerr << "\n!!! Errror in evaluating the DAE function "
                 "fcn(t,x,f,param).\n";
    std::cerr << "      " << exc.what() << '\n';
  }

  // extract the SigmaVector matrix from the SigmaVector vector f
  for (int i = 0; i < n_; i++) {
    // auto loop through key/value pairs (kv) from each f[i]
    for (auto& kv : f[i].sigmaMap()) {
      sigma_matrix_[i] = f[i].sigmaMap();
      nnz++;
    }
    if (i % 1000 == 0 && i > 0)
      printf("i = %i\n", i);
  }
  printf("Computed Sparse; number non-zero: %i\n", nnz);
  delete[] y;
  delete[] f;
}

void SigmaMatrix::printSparseSigmaMatrix() const {
  // Not the cleanest (or most intelligent) implementation of this, but this
  // should never really be called when N > ~10 so it doesn't matter.
  for (int i = 0; i < n_; i++) {  // rows
    for (int j = 0; j < n_; j++)  // columns
      if (sigma_matrix_[i].find(j) != sigma_matrix_[i].end())
        printf(" %i ", sigma_matrix_[i][j]);
      else
        printf(" - ");
    printf("\n");
  }
}