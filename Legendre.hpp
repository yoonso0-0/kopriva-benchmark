#include <array>
#include <cmath>
#include <assert.h>

#include "Constants.hpp"

#pragma once

void legendre_findroot_newton(const int polynomial_degree, const double x_init, double& x_root);
void legendre_pt_val_and_deriv(const int polynomial_degree, const double x_eval, double& f, double& df);

template <typename T, const size_t N>
void legendre_gauss_nodes_weights(std::array<T, N>& nodes, std::array<T, N>& weights)
{
/*
 *  Calculate gauss points and quadrature weights
 *  
 *  Parameters
 *  ----------
 *  Output, array<double, N> nodes : Gauss-Legendre quadrature points (degree : N)
 *  Output, array<double, N> weights : quadrature weights (degree : N)
 */
    assert( N >= 1 && "legendre_gauss_nodes_weights - input array has zero length");
    
    if (N==1)
    {
        nodes[0] = 0.0;
        weights[0] = 2.0;
        return;
    }
    else if (N==2)
    {
        nodes[0] = -1.0/sqrt(3.);
        nodes[1] = -nodes[0];
        weights[0] = 1.0;
        weights[1] = weights[0];
        return;
    }
    else
    {
        double delta_x;
        double Lx;
        double dLx;

        for (size_t j = 0; j < N/2; j++){

            // initial guess - Chebyshev Gauss points (Kopriva p.63)
            double x_j = - cos( M_PI * (2*j+1) / (2*N) );
        
            for(size_t i_root = 0; i_root <= constants::nmax_iter_rootsolving ; i_root++)
            {
                legendre_pt_val_and_deriv(N, x_j, Lx, dLx);

                delta_x = - Lx / dLx;
                x_j = x_j + delta_x;

                // std::cout << "iter = " << i_root << "\t x = " << x_j << "\t fx = " << Lx << std::endl;

                if ( std::abs(delta_x) < constants::epsilon_double )
                {
                    break;
                }
            }
            nodes[j] = x_j;
            nodes[N-1-j] = - nodes[j];
            weights[j] = 2.0 / ((1.0 - pow(x_j, 2.0)) * pow(dLx, 2.0));
            weights[N-1-j] = weights[j];
        }
        if (N % 2) {
            legendre_pt_val_and_deriv(N, 0., Lx, dLx);
            nodes[N/2] = 0.0;
            weights[N/2] = 2.0 / pow(dLx, 2.0);
        }
        return;
    }

    assert(0 && "You should not reach here!");
} // void legendre_gauss_nodes_weights

template <typename T, size_t N>
void legendre_gauss_lobatto_nodes_weights(std::array<T, N>, std::array<T, N>)
{
    // ...
}
