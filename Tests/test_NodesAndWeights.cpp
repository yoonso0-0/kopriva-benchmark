#include <iomanip>
#include <iostream>

#include <boost/numeric/ublas/vector.hpp>

#include "../NodesAndWeights.hpp"

namespace tests {

void test_gauss_nodes_weights() {
    //
    // Benchmark for N+1=7 (Kopriva p.67 Table 3.1)
    //
    const size_t benchmark_degree = 7;
    std::cout << " < Benchmark for N=" << benchmark_degree - 1 << " (Kopriva p.67 Table 3.1) > "
              << std::endl;

    boost::numeric::ublas::vector<double> nodes(benchmark_degree);
    boost::numeric::ublas::vector<double> weights(benchmark_degree);

    //
    // Gauss points
    //
    std::cout << " * Gauss" << std::endl;
    legendre::gauss_nodes_and_weights(nodes, weights);

    std::cout << std::setw(3) << "j" << std::setw(20) << "x_j" << std::setw(20) << "w_j"
              << std::endl;
    for (size_t i = 0; i < benchmark_degree / 2 + 1; i++) {
        std::cout << std::setw(3) << i << std::fixed << std::setprecision(15) << std::setw(20)
                  << nodes(i) << std::setw(20) << weights(i) << std::endl;
    }

    //
    // Gauss-Lobatto points
    //
    std::cout << " * Gauss-Lobatto" << std::endl;
    legendre::gauss_lobatto_nodes_and_weights(nodes, weights);

    std::cout << std::setw(3) << "j" << std::setw(20) << "x_j" << std::setw(20) << "w_j"
              << std::endl;
    for (size_t i = 0; i < benchmark_degree / 2 + 1; i++) {
        std::cout << std::setw(3) << i << std::fixed << std::setprecision(15) << std::setw(20)
                  << nodes(i) << std::setw(20) << weights(i) << std::endl;
    }
}

}  // namespace tests
