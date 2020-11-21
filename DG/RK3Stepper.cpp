#include <boost/numeric/ublas/vector.hpp>

#include <boost/numeric/ublas/io.hpp>

#include "../Constants.hpp"
#include "NodalDG.hpp"

#include "RK3Stepper.hpp"

void RK3Stepper(NodalDG &dg, const double time_step_size, const double t_n,
                double (*boundary_condition_function)(double)) {
    /*
     *  RK3 time integration
     *  (Kopriva Algorithm 62)
     */
    double t;
    boost::numeric::ublas::vector<double> G_j(dg.N_);
    std::fill(G_j.begin(), G_j.end(), 0.);

    for (size_t m = 0; m < 3; m++) {
        t = t_n + constants::RK3_bms[m] * time_step_size;
        dg.dg_timederivative(boundary_condition_function(t));

        for (size_t j = 0; j < dg.N_; j++) {
            G_j[j] = constants::RK3_ams[m] * G_j[j] + dg.solution_time_derivative_[j];
            dg.solution_array_[j] += constants::RK3_gms[m] * time_step_size * G_j[j];
        }
    }
    
    // Dirichlet boundary condition on the left
    dg.solution_array_[0] = boundary_condition_function(t_n + time_step_size);

}  // void RK3Stepper
