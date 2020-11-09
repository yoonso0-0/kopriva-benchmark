#include <assert.h>
#include <cmath>

#include <boost/numeric/ublas/vector.hpp>

#include "Constants.hpp"
#include "Legendre.hpp"

#pragma once

namespace legendre {

template <typename T>
void gauss_nodes_and_weights(boost::numeric::ublas::vector<T> &nodes,
                             boost::numeric::ublas::vector<T> &weights) {
  /*
   *  Calculate Gauss points and quadrature weights
   *  (Kopriva Algorithm 23)
   *
   *  Parameters
   *  ----------
   *  Output (implicit), vector<double> nodes (N) : Gauss quadrature points (degree : N)
   *  Output (implicit), vector<double> weights (N) : quadrature weights (degree : N)
   *
   */
  assert((nodes.size() == weights.size()) &&
         "legendre::gauss_nodes_weights - nodes and weights have different length");
  assert(nodes.size() >= 1 && "legendre::gauss_nodes_weights - input vector has zero length");

  size_t N = nodes.size();

  if (N == 1) {
    nodes[0] = 0.0;
    weights[0] = 2.0;
    return;
  } else if (N == 2) {
    nodes[0] = -1.0 / sqrt(3.);
    weights[0] = 1.0;
    nodes[1] = -nodes[0];
    weights[1] = weights[0];
    return;
  } else {
    T delta_x;
    T Lx;
    T dLx;

    for (size_t j = 0; j < N / 2; j++) {
      // initial guess - Chebyshev Gauss points
      T x_j = -cos(M_PI * (2 * j + 1) / (2 * N));

      // find root via Newton-Raphson method
      for (size_t i_root = 0; i_root <= constants::nmax_iter_rootsolving; i_root++) {
        legendre::polynomial_and_derivative(N, x_j, Lx, dLx);
        delta_x = -Lx / dLx;
        x_j = x_j + delta_x;

        if (std::abs(delta_x) < constants::epsilon_double)
          break;
      }
      legendre::polynomial_and_derivative(N, x_j, Lx, dLx);
      nodes[j] = x_j;
      nodes[N - 1 - j] = -nodes[j];
      weights[j] = 2.0 / ((1.0 - pow(x_j, 2.0)) * pow(dLx, 2.0));
      weights[N - 1 - j] = weights[j];
    }

    // if polynomial degree N is odd, manually assign weight for x=0 node
    if (N % 2) {
      legendre::polynomial_and_derivative(N, 0., Lx, dLx);
      nodes[N / 2] = 0.0;
      weights[N / 2] = 2.0 / pow(dLx, 2.0);
    }
    return;
  }
}  // void legendre_gauss_nodes_weights

template <typename T>
void gauss_lobatto_nodes_and_weights(boost::numeric::ublas::vector<T> &nodes,
                                     boost::numeric::ublas::vector<T> &weights) {
  /*
   *  Calculate Gauss-Lobatto points and quadrature weights
   *  (Kopriva Algorithm 25)
   *
   *  Parameters
   *  ----------
   *  Output (implicit), vector<double> nodes (N) : Gauss-Lobatto quadrature points (degree : N)
   *  Output (implicit), vector<double> weights (N) : quadrature weights (degree : N)
   *
   */

  assert((nodes.size() == weights.size()) &&
         "legendre::gauss_lobatto_nodes_weights - nodes and weights have different length");
  assert(nodes.size() >= 2 &&
         "legendre::gauss_lobatto_nodes_weights - input vector length must be >= 2");

  size_t N = nodes.size();

  if (N == 2) {
    nodes[0] = -1.0;
    nodes[1] = 1.0;
    weights[0] = 1.0;
    weights[1] = weights[0];
  } else {
    // endpoints x = -1.0, x = 1.0
    nodes[0] = -1.0;
    nodes[N - 1] = 1.0;
    weights[0] = 2.0 / (N * (N - 1));
    weights[N - 1] = weights[0];

    T delta_x;
    T q_x;
    T dq_x;
    T L_x;

    for (size_t j = 1; j < N / 2; j++) {
      // initial guess (Kopriva eq 3.7)
      T x_j = -cos((M_PI * (j + 0.25) - 3.0 / (8.0 * M_PI * (j + 0.25))) / (N - 1));

      // find root via Newton-Raphson method
      for (size_t i_root = 0; i_root <= constants::nmax_iter_rootsolving; i_root++) {
        legendre::q_and_L(N - 1, x_j, q_x, dq_x, L_x);
        delta_x = -q_x / dq_x;
        x_j = x_j + delta_x;

        if (std::abs(delta_x) < constants::epsilon_double)
          break;
      }
      legendre::q_and_L(N - 1, x_j, q_x, dq_x, L_x);
      nodes[j] = x_j;
      nodes[N - 1 - j] = -nodes[j];
      weights[j] = 2.0 / (N * (N - 1) * pow(L_x, 2.0));
      weights[N - 1 - j] = weights[j];
    }

    // if polynomial degree N is odd, manually assign weight for x=0 node
    if (N % 2) {
      legendre::q_and_L(N - 1, 0.0, q_x, dq_x, L_x);
      nodes[N / 2] = 0.0;
      weights[N / 2] = 2.0 / (N * (N - 1) * pow(L_x, 2.0));
    }
    return;
  }
}  // void gauss_lobatto_nodes_weights

}  // namespace legendre
