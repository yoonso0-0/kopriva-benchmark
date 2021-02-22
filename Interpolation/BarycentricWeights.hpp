#include <boost/numeric/ublas/vector.hpp>

#pragma once

namespace interpolation {

template <typename T>
void barycentric_weights(
    boost::numeric::ublas::vector<T>& interpolation_points,
    boost::numeric::ublas::vector<T>& barycentric_weights) {
  /*
   *  Calculate barycentric weights for Lagrange interpolation
   *  (Kopriva Algorithm 30)
   *
   *  Parameters
   *  ----------
   *  Input, vector<double> interpolation_points (N)
   *  Output (implicit), vector<double> barycentric_weights (N)
   *
   */
  size_t N{interpolation_points.size()};
  assert((N == barycentric_weights.size()) &&
         "interpolation::barycentric_weights - input vector lengths not match");

  for (size_t j = 0; j < N; j++) {
    barycentric_weights[j] = 1.0;
  }
  for (size_t j = 1; j < N; j++) {
    for (size_t k = 0; k < j; k++) {
      barycentric_weights[k] =
          barycentric_weights[k] *
          (interpolation_points[k] - interpolation_points[j]);
      barycentric_weights[j] =
          barycentric_weights[j] *
          (interpolation_points[j] - interpolation_points[k]);
    }
  }
  for (size_t j = 0; j < N; j++) {
    barycentric_weights[j] = 1.0 / barycentric_weights[j];
  }
  return;
}  // void barycentric_weights

}  // namespace interpolation
