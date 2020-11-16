#include <boost/numeric/ublas/vector.hpp>

#include "../Constants.hpp"

namespace interpolation {

void LagrangeInterpolatingPolynomials(double x_eval, boost::numeric::ublas::vector<double> &nodes,
                                      boost::numeric::ublas::vector<double> &weights,
                                      boost::numeric::ublas::vector<double> &interpolated_values);

}  // namespace interpolation
