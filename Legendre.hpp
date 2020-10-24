#include <array>
#include <cmath>

#include "Constants.hpp"

#pragma once

namespace legendre{

void polynomial_and_derivative(const int poly_degree_N, const double x_eval, double& L_x, double& dL_x);
void q_and_L(const int poly_degree_N, const double x_eval, double& q_x, double& dq_x, double& L_x);

} // namespace legendre
