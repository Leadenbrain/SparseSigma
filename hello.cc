#include <chrono>
#include "sparseSigMat.h"

/*
Our DAE Function: f(i) = -h^2*dx(i)/dt + D*(x(i+1) - 2 * x(i) + x(i-1))
  -h = spacing
  -D = diffison constant
Simulates molecular diffusion

Input: time, variables, functions, and parameters for DAE
Output: Function is read into the compSigmaMatrix() function
Algorithm:
  -Initialize variables needed
  -Preprocessor define our function to shift and allow casting as function for
SigmaMatrixFCN
  -Define DAE using molecular diffusion equation
*/
// Molecular Diffusion Simulation
template <typename T>
void compDAEmol(const T* x, T* f, void* p) {
  // Initialize variables
  int n = *(int*)p;                // number of molecules is the only param
  const double x0{0}, xn{10};      // space interval
  const double h{(xn - x0) / 10};  // spacing
  const double h2{h * h};
  const double D{0.96};   // diffusion constant
  const double u_b{0.1};  // Boundary condition

// Shift function/allow casting as a function in code
#define f(i) f[i - 1]
#define u(i) x[i - 1]

  // First equation in DAE
  f(1) = -h2 * Diff(u(1), 1) + D * (u(2) - 2 * u(1) + u_b);

  // Loop through all equations in DAE but first and last
  for (int i = 2; i <= n; i++)
    f(i) = -h2 * Diff(u(i), 1) + D * (u(i + 1) - 2 * u(i) + u(i - 1));

  // Last equation in DAE
  f(n + 1) = u(n + 1) - u(n - 1);

#undef u
#undef f
}

// Layne Watson Simulation
template <typename T>
void compDAELW(const T* x, T* f, void* param) {
  // The param argument conveys the problem size:
  int n = *((int*)param);
  int nplus1 = n + 1;
  T lambda = x[0];

  // Define f[0] = function S that specifies s to be arc length:
  f[0] = -1.0;
  for (int k = 0; k < nplus1; k++)
    f[0] += sqr(Diff(x[k], 1));

  T sum = 0.0;
  for (int k = 1; k <= n; k++)
    sum += x[k];

  // The algebraic equations are f[1] to f[n]
  for (int i = 1; i <= n; i++) {
    f[i] = x[i] - lambda * exp(cos(i * sum));
  }
}

// Function to allow us to quickly demonstrate code capability
void runSparseSimulation(int N, SigmaMatrixFcn DAE_FUN) {
  SigmaMatrix* sig_mat = new SigmaMatrix(N + 1, DAE_FUN);
  auto start = std::chrono::high_resolution_clock::now();
  sig_mat->compSparseSigmaMatrix(DAE_FUN, &N);
  auto stop = std::chrono::high_resolution_clock::now();
  auto duration =
      std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
  printf("Time to run compSparseSigma w/ %i mol: %f s\n", N,
         duration.count() / 1000000.0f);
  if (N <= 10)
    sig_mat->printSparseSigmaMatrix();

  delete sig_mat;
}

/*
Just a simple main function to run our program.
This code will ultimately be added to a C++ library. Instead of running
independently.
*/
int main(void) {
  // For sparse demonstration, we need over 100k (for our purposes)
  // 1 million is more than sufficient to ensure code is good
  int NSparse = 50000000;
  // For large N, we cant see the results easily. Therefore we will run the
  // molecular diffusion simulation at 10 to check results
  // ->results should be tridiagonal, 1 in diagonal, 0 in off diagonal
  // Then we will run another simulation that is dense in nature; the Layne
  // Watson Mechanism
  // -> Results should show 1 in every entry of first row, then zero in every
  // entry after
  int NDense = 10;
  runSparseSimulation(NSparse, compDAEmol);
  runSparseSimulation(NDense, compDAEmol);
  runSparseSimulation(NDense, compDAELW);
}
