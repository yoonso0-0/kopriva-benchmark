#include <iostream>
#include <iomanip>

#include "Constants.hpp"
#include "Legendre.hpp"
#include "NodesAndWeights.hpp"
#include "Interpolation/BarycentricWeights.hpp"

void benchmark_gauss_nodes_weights();
void test_legendre_pt_val_and_deriv();


// =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+= //
//                                     MAIN
// =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+= //
int main(int argc, char const *argv[])
{
    // test_legendre_pt_val_and_deriv();
    // benchmark_gauss_nodes_weights();

    std::array<double, 4> pts {-2, -1, 1, 2};
    std::array<double, 4> weights;

    interpolation::barycentric_weights(pts, weights);

    for (size_t i = 0; i < 4; i++)
    {
        std::cout << pts[i] << "   " << weights[i] << std::endl;
    }
    

    return 0;
}

// -------------------------------- TEST ROUTINES -------------------------------- //

void test_legendre_pt_val_and_deriv()
{
    std::cout << "* ------ TEST : P_N(x) ------ *" << std::endl;
    double f, df;
    legendre::polynomial_and_derivative(0, 0.4, f, df);

    std::cout << f << std::endl;
    std::cout << df << std::endl;
}

void benchmark_gauss_nodes_weights()
{
    // 
    // Benchmark for N+1=7 (Kopriva p.67 Table 3.1)
    // 
    const size_t benchmark_degree = 7;
    std::cout << " < Benchmark for N=" << benchmark_degree-1 << " (Kopriva p.67 Table 3.1) > " << std::endl;
    
    std::array<double, benchmark_degree> nodes;
    std::array<double, benchmark_degree> weights;

    // print out results up to 15-th precision
    std::cout << std::setprecision(15);
    std::cout << std::fixed;

    // 
    // Gauss points
    // 
    std::cout << " * Gauss" << std::endl;
    legendre::gauss_nodes_and_weights(nodes, weights);

    std::cout << std::setw(3) << "j" << std::setw(20) << "x_j" << std::setw(20) << "w_j" << std::endl;
    for (size_t i = 0; i < benchmark_degree/2 + 1; i++)
    {
        std::cout << std::setw(3) << i 
                  << std::setw(20) << nodes.at(i)
                  << std::setw(20) << weights.at(i)
                  << std::endl;
    }
    
    // 
    // Gauss-Lobatto points
    // 
    std::cout << " * Gauss-Lobatto" << std::endl;
    legendre::gauss_lobatto_nodes_and_weights(nodes, weights);
    
    std::cout << std::setw(3) << "j" << std::setw(20) << "x_j" << std::setw(20) << "w_j" << std::endl;
    for (size_t i = 0; i < benchmark_degree/2 + 1; i++)
    {
        std::cout << std::setw(3) << i 
                  << std::setw(20) << nodes.at(i)
                  << std::setw(20) << weights.at(i)
                  << std::endl;
    }

}
