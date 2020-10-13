#include <iostream>

#include "Constants.hpp"
#include "Legendre.hpp"

// 
void test_legendre_findroot_newton();
void test_legendre_pt_val_and_deriv();

// =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+= //
//                                     MAIN
// =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+= //
int main(int argc, char const *argv[])
{
    
    // test_legendre_pt_val_and_deriv();
    // test_legendre_findroot_newton();

    int a[3] {1,2,3};
    legendre_gauss_nodes_weights(a, a);

    return 0;
}

// -------------------------------- TEST ROUTINES -------------------------------- //
void test_legendre_findroot_newton()
{
    std::cout << "* ------ TEST : Root of Pn(x) ------ *" << std::endl;
    int stat;
    double root;
    legendre_findroot_newton(7, -0.004, root, stat);
    std::cout << "root = " << root << std::endl;
    std::cout << "stat = " << stat << std::endl;
}

void test_legendre_pt_val_and_deriv()
{
    std::cout << "* ------ TEST : P_N(x) ------ *" << std::endl;
    int stat;
    double f, df;
    legendre_pt_val_and_deriv(4, 0.4, f, df, stat);

    std::cout << f << std::endl;
    std::cout << df << std::endl;
}
