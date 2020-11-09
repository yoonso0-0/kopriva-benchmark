#include <cmath>
#include <cstddef>

#pragma once

namespace constants {

constexpr double epsilon_double{1e-16};
constexpr size_t nmax_iter_rootsolving{100};

inline int AlmostEqual(double a, double b) {
    double abs_diff = abs(a - b);

    if ((a == 0.) or (b == 0.)) {
        if (abs_diff < epsilon_double) {
            return 1;
        } else {
            return 0;
        }
    } else {
        if ((abs_diff < epsilon_double * abs(a)) and (abs_diff < epsilon_double * abs(b))) {
            return 1;
        } else {
            return 0;
        }
    }
}

}  // namespace constants
