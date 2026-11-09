#pragma once

#include <utility>
#include <vector>
#include <type_traits>
#include <stdexcept>
#include <concepts>
#include <fstream>

template<typename T>
concept arithmetic = std::is_arithmetic<T>::value;

namespace numerical {
template<arithmetic T>
class Matrix {
private:
    size_t n, m;
    std::vector<T> v;
    size_t get_ind(size_t i, size_t j) const { return i * m + j; }

public:
    Matrix() = default;
    Matrix(size_t n_, size_t m_) : n(n_), m(m_), v(n_ * m_, T(0)) {}

    void read_from_file(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open file: " + filename);
        }
        size_t rows, cols;
        if (!(file >> rows >> cols)) throw std::runtime_error("Cannot read matrix dimensions");
        Matrix temp(rows, cols);
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                if (!(file >> temp(i, j))) {
                    throw std::runtime_error("Cannot read matrix element");
                }
            }
        }
        *this = std::move(temp);
    }

    T& operator() (size_t i, size_t j) { return v[get_ind(i, j)]; }
    const T& operator() (size_t i, size_t j) const { return v[get_ind(i, j)]; }

    T& at(size_t i, size_t j) {
        if (i >= n || j >= m) throw std::out_of_range("Matrix index out of range");
        return v[get_ind(i, j)];
    }
    const T& at(size_t i, size_t j) const {
        if (i >= n || j >= m) throw std::out_of_range("Matrix index out of range");
        return v[get_ind(i, j)];
    }

    Matrix operator+(const Matrix& other) const {
        if (n != other.n || m != other.m) throw std::invalid_argument("Matrix dimensions mismatch");
        Matrix res(n, m);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < m; ++j) {
                res(i, j) = (*this)(i, j) + other(i, j);
            }
        }
        return res;
    }

    Matrix& operator+=(const Matrix& other) {
        if (n != other.n || m != other.m) throw std::invalid_argument("Matrix dimensions mismatch");
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < m; ++j) {
                (*this)(i, j) += other(i, j);
            }
        }
        return *this;
    }

    Matrix operator-(const Matrix& other) const {
        if (n != other.n || m != other.m) throw std::invalid_argument("Matrix dimensions mismatch");
        Matrix res(n, m);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < m; ++j) {
                res(i, j) = (*this)(i, j) - other(i, j);
            }
        }
        return res;
    }

    Matrix& operator-=(const Matrix& other) {
        if (n != other.n || m != other.m) throw std::invalid_argument("Matrix dimensions mismatch");
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < m; ++j) {
                (*this)(i, j) -= other(i, j);
            }
        }
        return *this;
    }

    Matrix operator*(const Matrix& other) const {
        if (m != other.n) throw std::invalid_argument("Matrix dimensions mismatch");
        Matrix res(n, other.m);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < other.m; ++j) {
                T sum = 0;
                for (size_t p = 0; p < m; ++p) {
                    sum += (*this)(i, p) * other(p, j);
                }
                res(i, j) = sum;
            }
        }
        return res;
    }

    Matrix& operator*=(const Matrix& other) { 
        *this = *this * other;
        return *this;
    }

    Matrix& operator*=(const T& val) {
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < m; ++j) {
                (*this)(i, j) *= val;
            }
        }
        return *this;
    }

    Matrix& operator/=(const T& val) {
        if (val == 0) throw std::runtime_error("Division by zero");
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < m; ++j) {
                (*this)(i, j) /= val;
            }
        }
        return *this;
    }

    Matrix operator*(const T& val) const {
        Matrix res = *this;
        res *= val;
        return res;
    }

    Matrix operator/(const T& val) const {
        Matrix res = *this;
        res /= val;
        return res;
    }

    size_t rows() const noexcept { return n; }
    size_t cols() const noexcept { return m; }
};
}