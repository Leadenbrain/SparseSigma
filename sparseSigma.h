// This class is going to be the workhorse of our algorithm
// We will run these sigmavectors through the DAE function, and using their
// overloaded operators will retrieve the order of each variable for a function.
// These will be combined together to make our sparse matrix
#ifndef INCLUDED_SPARSESIGMA_H
#define INCLUDED_SPARSESIGMA_H

#include <climits>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <limits>
#include <map>
#include <vector>

template <typename T>
using SparseVec = std::map<int, T>;
using SparseVecInt = SparseVec<int>;

class SigmaVector {
 public:
  //  Default empty
  SigmaVector() : val_(std::numeric_limits<double>::quiet_NaN()) {}
  // Overload with value
  explicit SigmaVector(int i) : val_(i) {}
  SigmaVector(double t) : val_(t) {}
  // Overload with value and map
  SigmaVector(const SigmaVector& t) : val_(t.val_), sigma_map_(t.sigma_map_) {}
  // Move operator to avoid temp value
  SigmaVector(SigmaVector&& t) : val_(t.val_) {
    sigma_map_ = std::move(t.sigma_map_);
  }
  SigmaVector& operator=(SigmaVector&& t) = default;

  // Allow us to retrieve the values at a key
  int operator[](unsigned int i) const {
    if (sigma_map_.find(i) != sigma_map_.end())
      return sigma_map_.at(i);
    else
      return 0;
  }

  // Binary operators
  friend SigmaVector operator+(const SigmaVector& t1, const SigmaVector& t2);
  friend SigmaVector operator-(const SigmaVector& t1, const SigmaVector& t2);
  friend SigmaVector operator*(const SigmaVector& t1, const SigmaVector& t2);
  friend SigmaVector operator/(const SigmaVector& t1, const SigmaVector& t2);

  SigmaVector& operator+=(const SigmaVector& t1);
  SigmaVector& operator-=(const SigmaVector& t1);
  SigmaVector& operator*=(const SigmaVector& t1);
  SigmaVector& operator/=(const SigmaVector& t1);

  friend SigmaVector pow(const SigmaVector& t1, double d);
  // TODO: Check this! Likely not making it into project, please ignore if you
  // see this
  friend SigmaVector pow(const SigmaVector& t1, const SigmaVector& t2);

  // unary operators
  friend SigmaVector operator+(const SigmaVector& t) { return t; }
  friend SigmaVector operator-(const SigmaVector& t) {
    if (t.is_map_const())
      return SigmaVector(-t.val_);
    return t;
  }

// All of these std operators have the same implementation
// A macro can be used to easily implement these operators with minimal coding
// This will replace OP with the std function we want to call
#define OPERATOR(OP)                            \
  friend SigmaVector OP(const SigmaVector& t) { \
    if (t.is_map_const())                       \
      return SigmaVector(std::OP(t.val_));      \
    return t;                                   \
  }

  // These are not needed for this project, but needed in general and so they
  // are included.
  OPERATOR(cos);
  OPERATOR(tan);
  OPERATOR(sqrt);
  OPERATOR(exp);
  OPERATOR(log);
  OPERATOR(asin);
  OPERATOR(acos);
  OPERATOR(atan);

  // This operator is not in std in a way that fits our macro
  // Will just write this out personally.
  friend SigmaVector sqr(const SigmaVector& t) {
    if (t.is_map_const())
      return SigmaVector(t.val_ * t.val_);
    return t;
  }

  friend SigmaVector Diff(const SigmaVector& t, int d);

  /// Print elements  s[0], s[1], s[size-1];
  friend std::ostream& operator<<(std::ostream& s, const SigmaVector& t);

  void initialize_vector() { sigma_map_.clear(); }

  void initialize_vector(int p) {
    initialize_vector();
    sigma_map_[p] = 0;
  }

  // Public func to access private variable val_
  double val() const { return val_; };
  SparseVecInt sigmaMap() const { return sigma_map_; }

 private:
  //  If size==0 then no variables, its a constant
  bool is_map_const() const { return sigma_map_.size() == 0; }
  // If constant and the value is 0, it must be zero.
  bool is_zero()
      const {  // cout << sigma_vec_.size() << "    " << val_ << endl;
    return (is_map_const() && fabs(val_) == 0);
  }

  // if a constants, store its value
  double val_;
  SparseVecInt sigma_map_;
};

// This is going to get used extensively, but needs to come after SigmaVector
// declaration, so I threw it here.
typedef void (*SigmaMatrixFcn)(const SigmaVector* y,
                               SigmaVector* f,
                               void* dae_params);

#endif
