#include <boost/numeric/ublas/vector.hpp>

#include "../Constants.hpp"
#include "LagrangeInterpolatingPolynomials.hpp"

namespace interpolation {

void LagrangeInterpolatingPolynomials(
    double x_eval, boost::numeric::ublas::vector<double>& nodes,
    boost::numeric::ublas::vector<double>& weights,
    boost::numeric::ublas::vector<double>& interpolated_values) {
  size_t N{interpolated_values.size()};
  assert((N == nodes.size()) &&
         "interpolation::barycentric_weights - input vector lengths not match");
  assert((N == weights.size()) &&
         "interpolation::barycentric_weights - input vector lengths not match");

  bool xMatchesNode{0};

  for (size_t j = 0; j < N; j++) {
    interpolated_values[j] = 0.0;
    if (constants::AlmostEqual(x_eval, nodes[j])) {
      interpolated_values[j] = 1.0;
      xMatchesNode = 1;
    }
  }

  if (xMatchesNode)
    return;

  double s{0};
  double t{};
  for (size_t j = 0; j < N; j++) {
    t = weights[j] / (x_eval - nodes[j]);
    interpolated_values[j] = t;
    s += t;
  }
  for (size_t j = 0; j < N; j++) {
    interpolated_values[j] = interpolated_values[j] / s;
  }
  return;

}  // void LagrangeInterpolatingPolynomials

}  // namespace interpolation
