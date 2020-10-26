#include <array>

#pragma once

namespace interpolation {

template <typename T, size_t N>
void barycentric_weights(const std::array<T, N>& interpolation_points, std::array<T, N>& barycentric_weights)
{
/* 
 *  Calculate barycentric weights for Lagrange interpolation
 *  (Kopriva Algorithm 30)
 *  
 *  Parameters
 *  ----------
 *  Input, array<double, N> interpolation_points
 *  Output (implicit), array<double, N> barycentric_weights
 * 
 */
    for (size_t j = 0; j < N; j++)
    {
        barycentric_weights[j] = 1.0;
    }
    for (size_t j = 1; j < N; j++)
    {
        for (size_t k = 0; k < j; k++)
        {
            barycentric_weights[k] = barycentric_weights[k] * (interpolation_points[k] - interpolation_points[j]);
            barycentric_weights[j] = barycentric_weights[j] * (interpolation_points[j] - interpolation_points[k]);
        }
    }
    for (size_t j = 0; j < N; j++)
    {
        barycentric_weights[j] = 1.0 / barycentric_weights[j];
    }
    return;
}

}
