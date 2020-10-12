#include <iostream>
#include <cmath>

#include "Legendre.hpp"
#include "Constants.hpp"

using namespace constants;

void legendre_pt_val_and_deriv(int N, double x, double& f, double& df, int& ierr)
{
/* 
 *  Evaluates P_N(x) and P'_N(x) for a given N and x 
 *  
 *  Parameters
 *  ----------
 *  Input, int N : order of polynomial
 * 
 *  Input, double X : point at which P_N, P'_N are evaluated
 * 
 *  Output, double f : P_N(x)
 * 
 *  Output, double df : P'_N(x)
 * 
 *  Output, int ierr : error flag
 *      0, no errors
 *      -1, x is out of domain [-1, 1]
 *      -2, N is negative
 * 
 */
    if (std::abs(x) > 1.0)
    {
        ierr = -1;
        std::cout << "legendre_pt_val_and_deriv - abs(X) > 1.0" << std::endl;
        return;
    }
    if (N==0)
    {
        f = 1.;
        df = 0.;
        return;
    }
    else if (N==1)
    {
        f = x;
        df = 1.;
        return;
    }
    else if (N<0)
    {
        ierr = -2;
        std::cout << "legendre_pt_val_and_deriv - N < 0" << std::endl;
        return;
    }
    else
    {
        double lkm2 {1.};
        double lkm1 {x};
        double dlkm2 {0.};
        double dlkm1 {1.};

        double lk {0.};
        double dlk {0.};

        for (int k = 2; k <= N; k++)
        {
            lk = ( lkm1 * x * (2 * k - 1) - lkm2 * (k - 1) ) / k ;
            dlk = dlkm2 + (2 * k - 1) * lkm1 ;
            
            // store current values for next step
            lkm2 = lkm1;
            lkm1 = lk;
            dlkm2 = dlkm1;
            dlkm1 = dlk;
        }
        ierr = 0;
        f = lk;
        df = dlk;
    }
}

void findroot_legendre_newton(int N, double xi, double& x, int& ierr)
{
/*
 *  Solve P_N(x0) = 0 for a given N via Newton's method
 *  
 *  Parameters
 *  ----------
 *  Input, int N : order of polynomial
 *  Input, double xi : initial guess 
 *  Output, double x : root
 *  Output, int* ierr : error flag
 *      0, no errors
 *      -1, root unobtainable
 * 
 */

    double delta_x {1.0};

    double f, df;
    int stat_leg;

    for(int i=0; i <= NMAX_ITER_ROOT ; i++)
    {
        legendre_pt_val_and_deriv(N, xi, f, df, stat_leg);
        if (stat_leg)
        {
            std::cout << "findroot_legendre_newton - x out of bound" << std::endl;
            ierr = -1;
            return;
        }

        delta_x = - f / df;
        x = xi + delta_x;
        xi = x;

        std::cout << "i = " << i << "\t x = " << x << "\t fx = " << f << std::endl;

        if ( std::abs(f) < epsilon_double )
        {
            ierr = 0;
            return;
        }
    }

    ierr = -1;
    return;

}
