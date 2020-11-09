#include <iomanip>
#include <iostream>

#include "Tests/tests.hpp"

int main() {
  tests::test_legendre_pt_val_and_deriv();
  std::cout << "" << std::endl;

  tests::test_gauss_nodes_weights();
  std::cout << "" << std::endl;

  tests::test_polynomial_derivative_matrix();
}
