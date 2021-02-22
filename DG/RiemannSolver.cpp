#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include "../Constants.hpp"

namespace ublas = boost::numeric::ublas;

void riemann_solver(const ublas::vector<double>& u_boundary_internal,
                    const ublas::vector<double>& u_boundary_external,
                    const double normal_1d,
                    ublas::vector<double>& numerical_flux) {
  /*
   *  Riemann solver for 1D scalar wave equation
   *
   *  Parameters
   *  ----------
   *  Input, vector<double> u_boundary_internal, u_boundary_external :
   *  Input, double normal_1d :
   *  Output (implicit), vector<double> numerical_flux :
   *
   *  Actions
   *  ----------
   *  Given two state vectors from each side of boundary interfaces, calculate
   *  numerical flux for 1D scalar wave equation
   *
   */

  double w_minus{};
  double w_plus{};

  w_minus = 0.5 * (normal_1d * u_boundary_internal[1] +
                   constants::wave_speed * u_boundary_internal[2]);
  w_plus = 0.5 * (normal_1d * u_boundary_external[1] -
                  constants::wave_speed * u_boundary_external[2]);

  numerical_flux[0] = 0.;
  numerical_flux[1] = normal_1d * constants::wave_speed * (w_minus - w_plus);
  numerical_flux[2] = w_minus + w_plus;
}