
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include "../NodesAndWeights.hpp"

#include "../Interpolation/BarycentricWeights.hpp"
#include "../Interpolation/LagrangeInterpolatingPolynomials.hpp"
#include "../Interpolation/PolynomialDerivativeMatrix.hpp"

#include "NodalDG.hpp"

namespace ublas = boost::numeric::ublas;

void NodalDG::construct() {
    /*
     *  Constructor for NodalDG class
     *  (Kopriva Algorithm 59)
     */

    // calculate quadrature weights and store values
    switch (quadrature_type_) {
        case 1:
            legendre::gauss_nodes_and_weights(nodes_, quadrature_weights_);
            break;
        case 2:
            legendre::gauss_lobatto_nodes_and_weights(nodes_, quadrature_weights_);
            break;
    }

    ublas::vector<double> barycentric_weights(N_);
    interpolation::barycentric_weights(nodes_, barycentric_weights);

    // store the values of Lagrange interpolation polynomials at both boundaries
    interpolation::LagrangeInterpolatingPolynomials(-1.0, nodes_, barycentric_weights,
                                                    basis_values_left_);
    interpolation::LagrangeInterpolatingPolynomials(+1.0, nodes_, barycentric_weights,
                                                    basis_values_right_);

    // Derivative matrix
    ublas::matrix<double> D_ij(N_, N_);
    interpolation::polynomial_derivative_matrix(nodes_, D_ij);
    for (size_t j = 0; j < N_; j++) {
        for (size_t i = 0; i < N_; i++) {
            derivative_matrix_hat_(i, j) =
                -D_ij(j, i) * quadrature_weights_(j) / quadrature_weights_(i);
        }
    }

}  // void NodalDG::construct

void NodalDG::dg_derivative(const double boundary_value_left, const double boundary_value_right) {
    /*
     *  Spatial derivative with DG approximation
     *  (Kopriva Algorithm 60)
     */

    solution_spatial_derivative_ = ublas::prod(derivative_matrix_hat_, solution_array_);

    for (size_t j = 0; j < N_; j++) {
        solution_spatial_derivative_(j) += (boundary_value_right * basis_values_right_(j) -
                                            boundary_value_left * basis_values_left_(j)) /
                                           quadrature_weights_(j);
    }

}  // void dg_derivative

void NodalDG::dg_timederivative(const double boundary_value_time) {
    /*
     *  Time derivative
     *  (Kopriva Algorithm 61)
     *  Note some differences since GLL points are currently used.
     */

    double boundary_value_left;
    double boundary_value_right;

    if (wave_speed_ >= 0.) {
        boundary_value_left = boundary_value_time;
        boundary_value_right = solution_array_[N_ - 1];
    } else {
        boundary_value_left = solution_array_[0];
        boundary_value_right = boundary_value_time;
    }

    dg_derivative(boundary_value_left, boundary_value_right);
    solution_time_derivative_ = -wave_speed_ * solution_spatial_derivative_;

}  // void dg_timederivative
