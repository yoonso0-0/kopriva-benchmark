#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include "BarycentricWeights.hpp"

namespace interpolation {

template <typename T>
void polynomial_derivative_matrix(boost::numeric::ublas::vector<T> &nodes,
                                  boost::numeric::ublas::matrix<T> &derivative_matrix)
/*
 *  Calculate differentiation matrix for a set of nodes
 *  (Kopriva Algorithm 37)
 *
 *  Parameters
 *  ----------
 *  Input, vector<double> nodes(N) : nodes
 *  Output (implicit), matrix<double> derivative_matrix(N,N) : differentiation matrix
 *
 */
{
  size_t n_nodes{derivative_matrix.size2()};  // number of columns in derivative_matrix

  assert((n_nodes == derivative_matrix.size1()) &&
         "polynomial_derivative_matrix - matrix is not square");
  assert((n_nodes == nodes.size()) && "polynomial_derivative_matrix - nrows != vector dimension");

  // calculate barycentric weights for D_ij computation
  boost::numeric::ublas::vector<T> weights(n_nodes);
  interpolation::barycentric_weights(nodes, weights);

  for (size_t i = 0; i < n_nodes; i++) {
    // set diagonal terms initially to be zero to use 'negative sum trick'
    derivative_matrix(i, i) = 0.0;

    for (size_t j = 0; j < n_nodes; j++) {
      if (j != i) {
        // Kopriva eq. (3.48)
        derivative_matrix(i, j) = weights(j) / (weights(i) * (nodes(i) - nodes(j)));
        // carrying out negative sum trick
        derivative_matrix(i, i) -= derivative_matrix(i, j);
      }
    }
  }

}  // void polynomial_derivative_matrix

}  // namespace interpolation
