#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include "../Constants.hpp"
#include "../Interpolation/BarycentricWeights.hpp"
#include "../Interpolation/LagrangeInterpolatingPolynomials.hpp"
#include "../Interpolation/PolynomialDerivativeMatrix.hpp"
#include "../NodesAndWeights.hpp"

#include "NodalDG.hpp"

namespace ublas = boost::numeric::ublas;

void NodalDG::construct() {
  /*
   *  Constructor for NodalDG class
   *  (Kopriva Algorithm 59)
   */

  // calculate nodes and quadrature weights
  switch (quadrature_type_) {
    case 1:
      legendre::gauss_nodes_and_weights(nodes_, quadrature_weights_);
      break;
    case 2:
      legendre::gauss_lobatto_nodes_and_weights(nodes_, quadrature_weights_);
      break;
  }

  // calculate barycentric weights for nodes, which is required for evaluating
  // Lagrange interpolating polynomials
  ublas::vector<double> barycentric_weights(N_);
  interpolation::barycentric_weights(nodes_, barycentric_weights);

  // store the values of Lagrange interpolation polynomials at boundaries
  interpolation::LagrangeInterpolatingPolynomials(
      -1.0, nodes_, barycentric_weights, basis_values_left_);
  interpolation::LagrangeInterpolatingPolynomials(
      +1.0, nodes_, barycentric_weights, basis_values_right_);

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

void NodalDG::dg_timederivative(
    const ublas::vector<double> numerical_flux_left,
    const ublas::vector<double> numerical_flux_right) {
  /*
   *  DG Time derivative for 1D scalar wave equation
   *
   *  Parameters
   *  ----------
   *  Input, vector<double> numerical_flux_left
   *  Input, vector<double> numerical_flux_right
   *
   *  Action
   *  ----------
   *  evaluate time derivative of state vector at each points with given
   * numerical fluxes at domain boundaries, and save that values to member
   * variable time_derivative_
   *
   */

  // evaluate flux F and source term S with state vector U
  for (size_t j = 0; j < N_; j++) {
    sol_F_(j, 0) = 0.;
    sol_F_(j, 1) = pow(constants::wave_speed, 2.0) * sol_U_(j, 2);
    sol_F_(j, 2) = sol_U_(j, 1);

    sol_S_(j, 0) = -sol_U_(j, 1);
    sol_S_(j, 1) = 0.;
    sol_S_(j, 2) = 0.;
  }

  // calculate volume term
  spatial_derivative_ = -ublas::prod(derivative_matrix_hat_, sol_F_);

  for (size_t m = 0; m < 3; m++) {
    for (size_t j = 0; j < N_; j++) {
      // add boundary terms
      spatial_derivative_(j, m) +=
          -(numerical_flux_right(m) * basis_values_right_(j) +
            numerical_flux_left(m) * basis_values_left_(j)) /
          quadrature_weights_(j);
      // add source terms and apply coordinate transform
      time_derivative_(j, m) =
          spatial_derivative_(j, m) / jacobian_1d_(j) + sol_S_(j, m);
    }
  }
}  // void NodalDG::dg_timederivative

void NodalDG::plugin_boundary_values() {
  /*
   *  Copy x=-1 and x=+1 values of sol_U_ to
   *  boundary_state_left_ and boundary_state_right_
   */
  for (size_t i = 0; i < 3; i++) {
    boundary_state_left_(i) = sol_U_(0, i);
    boundary_state_right_(i) = sol_U_(N_ - 1, i);
  }
}  // void NodalDG::plugin_boundary_values()
