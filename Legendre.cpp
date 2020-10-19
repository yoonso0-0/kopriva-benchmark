#include <iostream>
#include <cmath>
#include <assert.h>

#include "Legendre.hpp"
#include "Constants.hpp"

void legendre_pt_val_and_deriv(const int poly_degree_N, const double x_eval, double& L_x, double& dL_x)
{
/* 
 *  Evaluates L_N(x) and dL_N(x) for a given N and x 
 *  
 *  Parameters
 *  ----------
 *  Input, int poly_degree_N : order of polynomial
 *  Input, double x_eval : point at which L_N(x), dL_N(x) are evaluated
 *  Output (implicit), double L_x : L_N(x)
 *  Output (implicit), double dL_x : dL_N(x)
 * 
 */

    assert( abs(x_eval) <= 1.0 && "legendre_pt_val_and_deriv - abs(X) > 1.0");
    assert( poly_degree_N >= 0 && "legendre_pt_val_and_deriv - N < 0");

    if (poly_degree_N==0)
    {
        L_x = 1.;
        dL_x = 0.;
        return;
    }
    else if (poly_degree_N==1)
    {
        L_x = x_eval;
        dL_x = 1.;
        return;
    }
    else
    {
        double L_k_minus_2 {1.};
        double L_k_minus_1 {x_eval};
        double dL_k_minus_2 {0.};
        double dL_k_minus_1 {1.};

        double L_k {0.};
        double dL_k {0.};

        for (size_t k = 2; k <= poly_degree_N; k++)
        {
            L_k = ( L_k_minus_1 * x_eval * (2 * k - 1) - L_k_minus_2 * (k - 1) ) / k ;
            dL_k = dL_k_minus_2 + (2 * k - 1) * L_k_minus_1 ;
            
            // store current values for next step
            L_k_minus_2 = L_k_minus_1;
            L_k_minus_1 = L_k;
            dL_k_minus_2 = dL_k_minus_1;
            dL_k_minus_1 = dL_k;
        }
        L_x = L_k;
        dL_x = dL_k;
        return;
    }
}
