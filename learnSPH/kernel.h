#pragma once
#include <cassert>
#include <Eigen/Dense>

#include <CompactNSearch/CompactNSearch.h>

namespace learnSPH
{
    namespace Kernel
    {
        constexpr double PI = 3.14159265358979323846;

        double kernel(Eigen::Vector3d x_i, Eigen::Vector3d x_j, double h);
        Eigen::Vector3d kernelGrad(Eigen::Vector3d x_i, Eigen::Vector3d x_j, double h);
        
        double cubicSpline(const double q);
        double cubicSplineGrad(const double q);
    };
};
