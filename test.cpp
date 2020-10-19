#include <iostream>
#include <iomanip>

#include "Constants.hpp"
#include "Legendre.hpp"

void benchmark_gauss_nodes_weights();
void test_legendre_pt_val_and_deriv();


// =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+= //
//                                     MAIN
// =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+= //
int main(int argc, char const *argv[])
{
    // test_legendre_pt_val_and_deriv();
    // test_legendre_findroot_newton();

    benchmark_gauss_nodes_weights();
    

    return 0;
}

// -------------------------------- TEST ROUTINES -------------------------------- //
// void test_legendre_findroot_newton()
// {
//     std::cout << "* ------ TEST : Root of Pn(x) ------ *" << std::endl;
//     double root;
//     legendre_findroot_newton(7, -0.4, root);
//     std::cout << "root = " << root << std::endl;
// }

void test_legendre_pt_val_and_deriv()
{
    std::cout << "* ------ TEST : P_N(x) ------ *" << std::endl;
    double f, df;
    legendre_pt_val_and_deriv(0, 0.4, f, df);

    std::cout << f << std::endl;
    std::cout << df << std::endl;
}

void benchmark_gauss_nodes_weights()
{
    // 
    // Benchmark for N+1=7 (Kopriva p.67 Table 3.1)
    // 
    std::array<double, 7> nodes;
    std::array<double, 7> weights;
    legendre_gauss_nodes_weights(nodes, weights);

    // output formatting
    std::cout << std::setprecision(15);
    std::cout << std::fixed;
    std::cout << std::setw(3) << "j" << std::setw(20) << "x_j" << std::setw(20) << "w_j" << std::endl;

    for (size_t i=0; i<4; i++)
    {
        std::cout << std::setw(3) << i 
                  << std::setw(20) << nodes.at(i)
                  << std::setw(20) << weights.at(i)
                  << std::endl;
    }
}
