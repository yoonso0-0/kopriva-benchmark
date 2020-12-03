#include "Constants.hpp"
#include "DG/NodalDG.hpp"

void set_initial_condition(NodalDG &dg) {
    /*
     *  Set the state vector U at t=0 with the superposition of : 
     *  - Gaussian packet, sigma=0.2, x0 = -0.7, from left to right
     *  - Gaussian packet, sigma=0.1, x0 = -0.7, from right to left, 
     */

    for (size_t j = 0; j < dg.N_; j++) {

        double x0;
        double sigma;
        double xmx0;
        double psi;

        // packet 1
        x0 = -0.7;
        sigma = 0.2;
        xmx0 = dg.nodes_mapped_[j] - x0;
        psi = exp(-pow(xmx0 / sigma, 2.0));

        dg.sol_U_(j, 0) = psi;
        dg.sol_U_(j, 1) = -2.0 * psi * constants::wave_speed * xmx0 / (pow(sigma, 2.0));
        dg.sol_U_(j, 2) = -2.0 * psi * xmx0 / (pow(sigma, 2.0));

        // packet 2
        x0 = +0.7;
        sigma = 0.1;
        xmx0 = dg.nodes_mapped_[j] - x0;
        psi = 0.5 * exp(-pow(xmx0 / sigma, 2.0));   // half amplitude
        dg.sol_U_(j, 0) += psi;
        dg.sol_U_(j, 1) += 2.0 * psi * constants::wave_speed * xmx0 / (pow(sigma, 2.0));
        dg.sol_U_(j, 2) += -2.0 * psi * xmx0 / (pow(sigma, 2.0));
    }
}
