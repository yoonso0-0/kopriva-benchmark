
#pragma once

namespace tests {

void print_results_test_legendre_pt_val_and_deriv(int N, double px, double px_analytic, double dpdx,
                                                  double dpdx_analytic);

void test_gauss_nodes_weights();

void test_legendre_pt_val_and_deriv();

void test_polynomial_derivative_matrix();

}  // namespace tests
