#pragma once

#include "matrix.hpp"
#include <concepts>
#include <numeric>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace numerical {
template <std::floating_point T> class LU {
private:
  Matrix<T> lu;
  std::vector<size_t> perm;
  int sign = 1;
  T tolerance;

public:
  explicit LU(const Matrix<T> &a, T eps = 1e-12)
      : lu(a), perm(a.rows()), tolerance(eps) {
    if (a.rows() == 0 || a.rows() != a.cols())
      throw std::invalid_argument("LU requires a nonempty square matrix");
    std::iota(perm.begin(), perm.end(), 0);
    size_t n = a.rows();
    for (size_t k = 0; k < n; ++k) {
      size_t pivot = k;
      for (size_t i = k + 1; i < n; ++i) {
        if (std::abs(lu(i, k)) > std::abs(lu(pivot, k)))
          pivot = i;
      }
      if (std::abs(lu(pivot, k)) <= tolerance)
        throw std::runtime_error("Matrix is singular or nearly singular");
      if (pivot != k) {
        for (size_t j = 0; j < n; ++j)
          std::swap(lu(k, j), lu(pivot, j));
        std::swap(perm[k], perm[pivot]);
        sign = -sign;
      }
      for (size_t i = k + 1; i < n; ++i) {
        lu(i, k) /= lu(k, k);
        for (size_t j = k + 1; j < n; ++j)
          lu(i, j) -= lu(i, k) * lu(k, j);
      }
    }
  }

  Matrix<T> solve(const Matrix<T> &b) const {
    size_t n = lu.rows();
    if (b.rows() != n)
      throw std::invalid_argument("RHS row count mismatch");
    Matrix<T> x(n, b.cols());
    for (size_t col = 0; col < b.cols(); ++col) {
      for (size_t i = 0; i < n; ++i) {
        T sum = b(perm[i], col);
        for (size_t j = 0; j < i; ++j)
          sum -= lu(i, j) * x(j, col);
        x(i, col) = sum;
      }
      for (size_t i = n; i-- > 0;) {
        T sum = x(i, col);
        for (size_t j = i + 1; j < n; ++j)
          sum -= lu(i, j) * x(j, col);
        x(i, col) = sum / lu(i, i);
      }
    }
    return x;
  }

  T determinant() const {
    T det = static_cast<T>(sign);
    for (size_t i = 0; i < lu.rows(); ++i)
      det *= lu(i, i);
    return det;
  }

  Matrix<T> inverse() const {
    size_t n = lu.rows();
    Matrix<T> identity(n, n);
    for (size_t i = 0; i < n; ++i)
      identity(i, i) = T{1};
    return solve(identity);
  }
};
} // namespace numerical