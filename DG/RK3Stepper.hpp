#include "NodalDG.hpp"

#pragma once

void RK3Stepper(NodalDG &dg, const double time_step_size, const double time,
                double (*boundary_condition_function)(double));
