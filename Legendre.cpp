#include <assert.h>
#include <cmath>
#include <iostream>

#include "Constants.hpp"
#include "Legendre.hpp"

namespace legendre {

void polynomial_and_derivative(const int poly_degree_N, const double x_eval, double &L_x,
                               double &dL_x) {
    /*
     *  Evaluates Legendre polynomial and its derivagive for given N and x
     *  (Kopriva Algorithm 22)
     *
     *  Parameters
     *  ----------
     *  Input, int poly_degree_N : degree of Legendre polynomial
     *  Input, double x_eval : point at which L_N(x), dL_N(x) are evaluated
     *  Output (implicit), double L_x : L_N(x)
     *  Output (implicit), double dL_x : dL_N(x)
     *
     */

    assert(abs(x_eval) <= 1.0 && "legendre_pt_val_and_deriv - abs(X) > 1.0");
    assert(poly_degree_N >= 0 && "legendre_pt_val_and_deriv - N < 0");

    if (poly_degree_N == 0) {
        L_x = 1.;
        dL_x = 0.;
        return;
    } else if (poly_degree_N == 1) {
        L_x = x_eval;
        dL_x = 1.;
        return;
    } else {
        double L_k_minus_2{1.};
        double L_k_minus_1{x_eval};
        double dL_k_minus_2{0.};
        double dL_k_minus_1{1.};

        double L_k{};
        double dL_k{};

        for (size_t k = 2; k <= poly_degree_N; k++) {
            L_k = (L_k_minus_1 * x_eval * (2 * k - 1) - L_k_minus_2 * (k - 1)) / k;
            dL_k = dL_k_minus_2 + (2 * k - 1) * L_k_minus_1;

            // store current values for next step
            L_k_minus_2 = L_k_minus_1;
            L_k_minus_1 = L_k;
            dL_k_minus_2 = dL_k_minus_1;
            dL_k_minus_1 = dL_k;
        }
        L_x = L_k;
        dL_x = dL_k;
        return;
    }
}  // void polynomial_and_derivative

void q_and_L(const int poly_degree_N, const double x_eval, double &q_x, double &dq_x, double &L_x) {
    /*
     *  Evaluates q, q', L_N(x)
     *  (Kopriva Algorithm 24)
     *
     *  Parameters
     *  ----------
     *  Input, int poly_degree_N : degree of Legendre polynomial
     *  Input, double x_eval : point at which q_x, dq_x, L_x are evaluated
     *  Output (implicit), double q_x : q(x) = L_{N+1}(x) - L_{N-1}(x)
     *  Output (implicit), double dq_x : dq(x) = dL_{N+1}(x) - dL_{N-1}(x)
     *  Output (implicit), double L_x : L_N(x)
     *
     */

    assert(abs(x_eval) <= 1.0 && "legendre_q_and_L - abs(X) > 1.0");
    assert(poly_degree_N >= 0 && "legendre_q_and_L - N < 0");

    if (poly_degree_N == 0) {
        L_x = 1.;
        q_x = 0.;
        dq_x = 0.;
        return;
    } else if (poly_degree_N == 1) {
        L_x = x_eval;
        q_x = 1.;
        dq_x = 0.;
        return;
    } else {
        double L_k_minus_2{1.};
        double L_k_minus_1{x_eval};
        double dL_k_minus_2{0.};
        double dL_k_minus_1{1.};

        double L_k{};
        double dL_k{};

        for (size_t k = 2; k <= poly_degree_N; k++) {
            L_k = (L_k_minus_1 * x_eval * (2 * k - 1) - L_k_minus_2 * (k - 1)) / k;
            dL_k = dL_k_minus_2 + (2 * k - 1) * L_k_minus_1;

            // store current values for next step
            L_k_minus_2 = L_k_minus_1;
            L_k_minus_1 = L_k;
            dL_k_minus_2 = dL_k_minus_1;
            dL_k_minus_1 = dL_k;
        }

        // note that the value of L_{N-1}(x) is stored in L_k_minus_2 at this point
        L_x = L_k;
        q_x = (L_k * x_eval - L_k_minus_2) * (2 * poly_degree_N + 1) / (poly_degree_N + 1);
        dq_x = (2 * poly_degree_N + 1) * L_k;
        return;
    }
}  // void q_and_L

}  // namespace legendre
