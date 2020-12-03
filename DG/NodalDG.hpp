#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include "../Constants.hpp"

namespace ublas = boost::numeric::ublas;

#pragma once

class NodalDG {
    /*
     *  definition for Nodal DG class
     *  for solving 1D scalar wave equation
     *  ----------
     *   
     *  ...
     *  
     */
   public:
    // Expansion order; specified by input
    size_t N_;

    // Quadrature type (nodes)
    // 1 : GL (not implemented yet : @201116)
    // 2 : GLL
    int quadrature_type_{2};

    // quadrature nodes and weights
    ublas::vector<double> nodes_;
    ublas::vector<double> quadrature_weights_;

    // values of Lagrange interpolating polynomials at right / left boundary
    ublas::vector<double> basis_values_left_;
    ublas::vector<double> basis_values_right_;

    // state vector, flux, source
    ublas::matrix<double> sol_U_;
    ublas::matrix<double> sol_F_;
    ublas::matrix<double> sol_S_;

    // temporal and spatial derivative of sol_U_
    ublas::matrix<double> time_derivative_;
    ublas::matrix<double> spatial_derivative_;

    // state vector at the left and right boundary i.e. x=+1 and x=-1.
    ublas::vector<double> boundary_state_left_;
    ublas::vector<double> boundary_state_right_;

    // Derivative matrix \hat{D}_{ij}
    ublas::matrix<double> derivative_matrix_hat_;

    // 1D Jacobian for coordinate map
    ublas::vector<double> jacobian_1d_;

    // mapped coordinates of nodes
    ublas::vector<double> nodes_mapped_;

    NodalDG(const size_t input_N = 8) {
        N_ = input_N;

        nodes_.resize(N_);
        quadrature_weights_.resize(N_);
        basis_values_left_.resize(N_);
        basis_values_right_.resize(N_);

        sol_U_.resize(N_, 3);
        sol_F_.resize(N_, 3);
        sol_S_.resize(N_, 3);

        boundary_state_left_.resize(3);
        boundary_state_right_.resize(3);

        time_derivative_.resize(N_, 3);
        spatial_derivative_.resize(N_, 3);

        derivative_matrix_hat_.resize(N_, N_);
        jacobian_1d_.resize(N_);
        nodes_mapped_.resize(N_);
    }

    void construct();
    void dg_timederivative(const ublas::vector<double> numerical_flux_left,
                           const ublas::vector<double> numerical_flux_right);
    void plugin_boundary_values();

};  // class NodalDG
