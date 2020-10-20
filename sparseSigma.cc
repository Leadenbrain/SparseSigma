// This file implements the SigmaVector class. There is a lot going on in here,
// but none of it is too complicated. I'm hoping that this early intro will
// suffice for the sparsity of comments throughout the file.\
// To start, I define a function to determine the max map values given two maps
//    -> This will be used to overload our operators to return the proper values
// Next I define a macro for the basic arithmetic operators and overload them
// Following this, I overload the binary operators for our class
//    -> This is done to handle the various cases that may arise given the maps
// Finally, I define our power operator and Differential operator
#include "sparseSigma.h"
#include <cassert>

inline void maxMap(const SparseVecInt& u,
                   const SparseVecInt& v,
                   SparseVecInt& w) {
  for (auto& ukv : u) {
    if (v.find(ukv.first) != v.end()) {
      if (ukv.second > v.at(ukv.first)) {
        w[ukv.first] = ukv.second;
      } else {
        w[ukv.first] = v.at(ukv.first);
      };
    } else {
      w[ukv.first] = ukv.second;
    }
  }
}

// Define our new sparse operations using overloading:
#define SIGMA_OP(OP)                                                      \
  SigmaVector operator OP(const SigmaVector& t1, const SigmaVector& t2) { \
    SigmaVector temp(t1);                                                 \
    temp OP## = t2;                                                       \
    return temp;                                                          \
  }

SIGMA_OP(+)
SIGMA_OP(-)
SIGMA_OP(*)
SIGMA_OP(/)

// Sparse Multiplication Operator
SigmaVector& SigmaVector::operator*=(const SigmaVector& t) {
  // both constants
  if (this->is_map_const() && t.is_map_const()) {
    // both are constants
    assert(this->sigma_map_.size() == 0);
    assert(t.sigma_map_.size() == 0);
    // Just multiply the values that exist
    this->val_ *= t.val_;
    return *this;
  }

  // one is zero
  if (this->is_zero() || t.is_zero()) {
    // one of the operands is zero.
    assert(this->sigma_map_.size() == 0 || t.sigma_map_.size() == 0);
    // the result must be 0
    this->val_ = 0;
    // the vector for this may not be of size 0 so resize
    this->sigma_map_.clear();
    return *this;
  }

  if (this->is_map_const())  // t is not constant
  {
    assert(this->sigma_map_.size() == 0);
    assert(t.sigma_map_.size() != 0);
    // // the value of the constant does not matter. Just copy the vectors
    this->sigma_map_ = t.sigma_map_;
    return *this;
  }

  if (t.is_map_const())  // *this is not constant
  {
    assert(this->sigma_map_.size() != 0);
    assert(t.sigma_map_.size() == 0);
    return *this;
  }

  // now both are non-constant
  maxMap(t.sigma_map_, this->sigma_map_, this->sigma_map_);
  return *this;
}

// Sparse Division Operator
SigmaVector& SigmaVector::operator/=(const SigmaVector& t) {
  // both constants
  if (this->is_map_const() && t.is_map_const()) {
    // both are constants
    assert(this->sigma_map_.size() == 0);
    assert(t.sigma_map_.size() == 0);
    this->val_ /= t.val_;
    return *this;
  }

  // one is zero
  if (this->is_zero()) {
    assert(this->sigma_map_.size() == 0);
    return *this;
  }

  if (this->is_map_const())  // t is not constant
  {
    assert(this->sigma_map_.size() == 0);
    assert(t.sigma_map_.size() != 0);
    // the value of the constant does not matter. Just copy the vectors
    this->sigma_map_ = t.sigma_map_;
    return *this;
  }

  if (t.is_map_const())  // this is not constant
  {
    assert(this->sigma_map_.size() != 0);
    assert(t.sigma_map_.size() == 0);
    return *this;
  }

  // now both are non-constant
  maxMap(t.sigma_map_, this->sigma_map_, this->sigma_map_);
  return *this;
}

// Sparse Addition Operator
SigmaVector& SigmaVector::operator+=(const SigmaVector& t) {
  if (this->is_map_const() && t.is_map_const()) {
    // both are constants
    assert(this->sigma_map_.size() == 0);
    assert(t.sigma_map_.size() == 0);
    this->val_ += t.val_;
    return *this;
  }

  if (this->is_map_const())  // t is not constant
  {
    assert(this->sigma_map_.size() == 0);
    assert(t.sigma_map_.size() != 0);
    // the value of the constant does not matter. Just copy the vectors
    this->sigma_map_ = t.sigma_map_;
    return *this;
  }

  if (t.is_map_const())  // this is not constant
  {
    assert(this->sigma_map_.size() != 0);
    assert(t.sigma_map_.size() == 0);
    assert(this->sigma_map_.size() != 0);
    assert(t.sigma_map_.size() == 0);
    return *this;
  }

  // now both are non-constant
  maxMap(t.sigma_map_, this->sigma_map_, this->sigma_map_);
  return *this;
}

// Sparse Subtraction Operator
SigmaVector& SigmaVector::operator-=(const SigmaVector& t) {
  if (this->is_map_const() && t.is_map_const()) {
    // both are constants
    assert(this->sigma_map_.size() == 0);
    assert(t.sigma_map_.size() == 0);
    this->val_ -= t.val_;
    return *this;
  }

  if (this->is_map_const())  // t is not constant
  {
    assert(this->sigma_map_.size() == 0);
    assert(t.sigma_map_.size() != 0);
    // the value of the constant does not matter. Just copy the vectors
    this->sigma_map_ = t.sigma_map_;
    return *this;
  }

  if (t.is_map_const())  // *this is not constant
  {
    assert(this->sigma_map_.size() != 0);
    assert(t.sigma_map_.size() == 0);
    return *this;
  }

  // now both are non-constant
  maxMap(t.sigma_map_, this->sigma_map_, this->sigma_map_);
  return *this;
}

// Sparse Power operation for single valued operator
SigmaVector pow(const SigmaVector& t1, double d) {
  // t1^0 = 1
  if (d == 0)
    return SigmaVector(1);

  // constant ^ d
  if (t1.is_map_const())
    return SigmaVector(std::pow(t1.val_, d));

  return t1;
}
// Sparse Power Operator of two vectors
// TODO: Check this
SigmaVector pow(const SigmaVector& t1, const SigmaVector& t2) {
  // t1^0
  if (t2.is_zero())
    return SigmaVector(1);

  // t1^constant
  if (t2.is_map_const())
    return t1;

  // constant^t2
  if (t1.is_map_const())
    return t2;

  assert(!t1.is_map_const() && !t2.is_map_const());
  SigmaVector tmp(t1);
  maxMap(t1.sigma_map_, t2.sigma_map_, tmp.sigma_map_);
  return tmp;
}

//  Differential operator, it will return presence of variable, plus its order
SigmaVector Diff(const SigmaVector& t, int d) {
  if (t.is_map_const()) {
    return SigmaVector(0);
  }

  // TODO: Check if I can do this with a std::move operator, might save the tmp
  // data, already hitting huge matrices though, so Im not disappointed
  SigmaVector tmp(t);
  assert(t.sigma_map_.size() != 0);
  for (auto& kv : t.sigma_map_)
    tmp.sigma_map_[kv.first] = kv.second + d;

  return tmp;
}
