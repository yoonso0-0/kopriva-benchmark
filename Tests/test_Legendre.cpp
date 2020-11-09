#include <cmath>
#include <iomanip>
#include <iostream>

#include "../Legendre.hpp"
#include "tests.hpp"

namespace tests {

void print_results_test_legendre_pt_val_and_deriv(int N, double x, double px, double px_analytic,
                                                  double dpdx, double dpdx_analytic) {
    std::cout << "N = " << N << ") "
              << "x = " << std::fixed << std::setprecision(2) << x << std::endl;
    std::cout << std::setw(10) << "PN   " << std::fixed << std::setprecision(15) << std::setw(20)
              << px_analytic << std::setw(20) << px << std::scientific << std::setprecision(4)
              << std::setw(15) << px / px_analytic - 1.0 << std::endl;

    std::cout << std::setw(10) << "PN'  " << std::fixed << std::setprecision(15) << std::setw(20)
              << dpdx_analytic << std::setw(20) << dpdx << std::scientific << std::setprecision(4)
              << std::setw(15) << dpdx / dpdx_analytic - 1.0 << std::endl;
}

void test_legendre_pt_val_and_deriv() {
    std::cout << " < test - P_N and P_N' > " << std::endl;

    std::cout << std::setw(7) << "" << std::setw(20) << "Analytic" << std::setw(20) << "Numerical"
              << std::setw(15) << "Err" << std::endl;

    double N, x, px, dpdx;
    double px_analytic, dpdx_analytic;

    // N = 0 should yield constant
    N = 0;
    x = 0.4;
    px_analytic = 1.0;
    dpdx_analytic = 0.0;
    legendre::polynomial_and_derivative(N, x, px, dpdx);
    print_results_test_legendre_pt_val_and_deriv(N, x, px, px_analytic, dpdx, dpdx_analytic);

    // N = 1 : x
    N = 1;
    x = 0.3;
    px_analytic = x;
    dpdx_analytic = 1.0;
    legendre::polynomial_and_derivative(N, x, px, dpdx);
    print_results_test_legendre_pt_val_and_deriv(N, x, px, px_analytic, dpdx, dpdx_analytic);

    // N = 3 : (5x^3 - 3x) / 2
    N = 3;
    x = 0.7;
    px_analytic = 2.5 * pow(x, 3.0) - 1.5 * x;
    dpdx_analytic = 7.5 * pow(x, 2.0) - 1.5;
    legendre::polynomial_and_derivative(N, x, px, dpdx);
    print_results_test_legendre_pt_val_and_deriv(N, x, px, px_analytic, dpdx, dpdx_analytic);

    // N = 4 : (35x^4 - 30x^2 + 3) / 8
    // ....
}

}  // namespace tests
