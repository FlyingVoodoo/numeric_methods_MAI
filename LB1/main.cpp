#include <iostream>
#include <fstream>
#include <iomanip>
#include <utility>
#include <stdexcept>
#include <string>

#include "include/lu.hpp"
#include "include/matrix.hpp"

constexpr double EPS = 1e-10;
constexpr size_t MAX_ITERATIONS = 100000;

template <typename T> void print_matrix(const numerical::Matrix<T>& A) {
  for (size_t i = 0; i < A.rows(); ++i) {
    for (size_t j = 0; j < A.cols(); ++j) {
      std::cout << std::setw(15) << A(i, j) << ' ';
    }
    std::cout << '\n';
  }
}

[[nodiscard]] std::pair<numerical::Matrix<double>, numerical::Matrix<double>> read_system(const std::string& filename) {
  std::ifstream file(filename);
  if (!file) throw std::runtime_error("Cannot open file: " + filename);
  size_t n, m;
  if (!(file >> n >> m)) throw std::runtime_error("Cannot read matrix dimensions");
  if (n == 0 || n != m) throw std::invalid_argument("Matrix must be square and n must be greater than 0");
  numerical::Matrix<double> A(n, m);
  numerical::Matrix<double> b(n, 1);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      if (!(file >> A(i, j))) throw std::runtime_error("Cannot read matrix");
    }
  }
  for (size_t i = 0; i < n; ++i) {
    if (!(file >> b(i, 0))) throw std::runtime_error("Cannot read RHS");
  }
  return {A, b};
}

void run(const std::string& name, auto&& func) {
  std::cout << "\n========== " << name << " ==========\n";
  try {
    func();
  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << '\n';
  }
}

auto main() -> signed {
  std::cout << std::fixed << std::setprecision(10);

  run("LU decomposition", [] {
    auto [A, b] = read_system("data/lu.txt");
    numerical::LU<double> lu(A);
    numerical::Matrix<double> x = lu.solve(b);
    std::cout << "Solution:\n";
    print_matrix(x);
    std::cout << "\nDeterminant: " << lu.determinant() << '\n';
    std::cout << "\nInverse matrix:\n";
    print_matrix(lu.inverse());
  });

  return 0;
}