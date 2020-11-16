
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

namespace ublas = boost::numeric::ublas;

#pragma once

class NodalDG {
    /*
     *   DG Class definition (Kopriva algorithm 58)
     */
   public:
    size_t N_;

    // Quadrature type (nodes)
    // 1 : GL (not implemented yet : @201116)
    // 2 : GLL
    int quadrature_type_{2};

    // wave speed
    double wave_speed_{1.0};

    // quadrature nodes and weights
    ublas::vector<double> nodes_;
    ublas::vector<double> quadrature_weights_;

    // values of Lagrange interpolating polynomials at right / left boundary
    ublas::vector<double> basis_values_left_;
    ublas::vector<double> basis_values_right_;

    // solution vector and its derivativs (Phi in Kopriva book)
    ublas::vector<double> solution_array_;
    ublas::vector<double> solution_spatial_derivative_;
    ublas::vector<double> solution_time_derivative_;

    // Derivative matrix \hat{D}_{ij}
    ublas::matrix<double> derivative_matrix_hat_;

    NodalDG(const size_t input_N = 8) {
        N_ = input_N;

        // Set vector lengths with specified expansion order N, given as an input
        nodes_.resize(N_);
        quadrature_weights_.resize(N_);
        basis_values_left_.resize(N_);
        basis_values_right_.resize(N_);
        solution_array_.resize(N_);
        solution_spatial_derivative_.resize(N_);
        solution_time_derivative_.resize(N_);
        derivative_matrix_hat_.resize(N_, N_);
    }

    void construct();
    void dg_derivative(const double boundary_value_left, const double boundary_value_right);
    void dg_timederivative(const double boundary_value_time);

};  // class NodalDG
