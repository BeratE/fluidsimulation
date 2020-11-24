#pragma once
#include <cassert>
#include <Eigen/Dense>

#include <CompactNSearch/CompactNSearch.h>

namespace learnSPH {
    namespace Kernel {
        namespace Parameter {
            constexpr double TUNING = 1.2;
            constexpr double EPSILON = 0.5;
        } 

        // Cubic Spline Kernel
        namespace CubicSpline {
            double support(double h);
            double weight(Eigen::Vector3d x_i, Eigen::Vector3d x_j, double h);
            Eigen::Vector3d gradWeight(Eigen::Vector3d x_i, Eigen::Vector3d x_j, double h);
            double cubicSpline(const double q);
            double gradCubicSpline(const double q);
        } // namespace CubicSpline
    } // namespace Kernel
} // namespace learnSPH
