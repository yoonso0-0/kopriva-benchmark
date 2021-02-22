#include <boost/numeric/ublas/vector.hpp>

#pragma once

void riemann_solver(const ublas::vector<double>& u_boundary_internal,
                    const ublas::vector<double>& u_boundary_external,
                    const double normal_1d,
                    ublas::vector<double>& numerical_flux);
