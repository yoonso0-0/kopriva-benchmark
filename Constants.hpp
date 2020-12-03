#include <cmath>
#include <cstddef>

#pragma once

namespace constants {

constexpr double epsilon_double{1e-16};
constexpr size_t nmax_iter_rootsolving{100};

constexpr double wave_speed{1.0};

// SSP-RK3 coefficients - Williamson's (Kopriva Table 4.1)
constexpr double RK3_ams[3] = {0.0, -5.0 / 9, -153.0 / 128};
constexpr double RK3_bms[3] = {0, 1.0 / 3, 3.0 / 4};
constexpr double RK3_gms[3] = {1.0 / 3, 15.0 / 16, 8.0 / 15};

inline int AlmostEqual(double a, double b) {
    double abs_diff = abs(a - b);

    if ((a == 0.) or (b == 0.)) {
        if (abs_diff < epsilon_double) {
            return 1;
        } else {
            return 0;
        }
    } else {
        if ((abs_diff < epsilon_double * abs(a)) and
            (abs_diff < epsilon_double * abs(b))) {
            return 1;
        } else {
            return 0;
        }
    }
}

}  // namespace constants
