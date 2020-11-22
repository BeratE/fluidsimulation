#include "kernel.h"
#define _USE_MATH_DEFINES
#include <math.h>
using namespace learnSPH::Kernel;

double CubicSpline::support(double h)
{
    return 2*h;
}

double CubicSpline::weight(Eigen::Vector3d x_i, Eigen::Vector3d x_j, double h)
{
    double q = (x_j - x_i).norm() / h;
    return pow(h, -3) * cubicSpline(q);
}

Eigen::Vector3d CubicSpline::gradWeight(Eigen::Vector3d x_i, Eigen::Vector3d x_j, double h)
{
    Eigen::Vector3d dist = (x_j - x_i);
    double q = dist.norm() / h;
    return pow(h, -4) * gradCubicSpline(q) * dist.normalized();
}

double CubicSpline::cubicSpline(const double q)
{
    assert(q >= 0.0);
        
    constexpr double alpha = 3.0 / (2.0 * M_PI);
    double value = 0.0;
        
    if (q < 1.0) {
        value = (2.0/3.0 - q*q + 0.5*q*q*q);
    }
    else if (q < 2.0) {
        value = pow((2.0-q), 3) / 6.0;
    }
        
    return alpha * value;
}

double CubicSpline::gradCubicSpline(const double q)
{
    assert(q >= 0.0);
        
    constexpr double alpha = 3.0 / (2.0 * M_PI);
    double value = 0.0;
        
    if (q < 1.0) {
        value = (3.0/2.0)*q*q - 2*q;
    }
    else if (q < 2.0) {
        value = -0.5*(2-q)*(2-q);
    }
        
    return alpha * value;
}

