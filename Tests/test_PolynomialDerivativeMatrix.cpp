#include <cstddef>
#include <iomanip>
#include <iostream>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include "../Interpolation/PolynomialDerivativeMatrix.hpp"

namespace tests {

void test_polynomial_derivative_matrix() {
    //
    // compare polynomial derivatives for a simple case:
    //  P(x)  = x^3 + 2x^2 + 5x + 1  (degree : 3)
    //  dP(x) = 3x^2 + 4x + 5
    //

    const size_t N = 4;

    boost::numeric::ublas::vector<double> test_points(N);
    test_points[0] = -0.5;
    test_points[1] = 0.0;
    test_points[2] = 0.31;
    test_points[3] = 0.4;
    // test_points[4] = 0.6;

    double dpdx_analytic[N];
    boost::numeric::ublas::vector<double> px_values(N);

    for (size_t i = 0; i < N; i++) {
        double x = test_points[i];
        //  P(x)  = x^3 + 2x^2 + 5x + 1
        px_values[i] = pow(x, 3.0) + 2.0 * pow(x, 2.0) + 5.0 * x + 1.0;
        //  dP(x) = 3x^2 + 4x + 5
        dpdx_analytic[i] = 3.0 * pow(x, 2.0) + 4.0 * x + 5.0;
    }

    boost::numeric::ublas::matrix<double> diff_matrix(N, N);
    interpolation::polynomial_derivative_matrix(test_points, diff_matrix);

    // store numerical derivatives
    boost::numeric::ublas::vector<double> dpdx_numerical(N);
    dpdx_numerical = boost::numeric::ublas::prod(diff_matrix, px_values);

    std::cout << " < test - derivative matrix > " << std::endl;

    std::cout << "Degree (N) = " << N << std::endl;

    std::cout << std::setw(20) << "Analytic" << std::setw(20) << "Numerical" << std::setw(15)
              << "Err" << std::endl;
    for (size_t i = 0; i < N; i++) {
        std::cout << std::fixed << std::setprecision(15) << std::setw(20) << dpdx_analytic[i]
                  << std::setw(20) << dpdx_numerical[i] << std::scientific << std::setprecision(4)
                  << std::setw(15) << dpdx_numerical[i] / dpdx_analytic[i] - 1.0 << std::endl;
    }
}

}  // namespace tests