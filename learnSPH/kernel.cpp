#include "kernel.h"

double learnSPH::Kernel::kernel(Eigen::Vector3d x_i, Eigen::Vector3d x_j,
                                double h)
{
    double q = (x_j - x_i).norm() / h;

    return pow(h, -3) * cubicSpline(q);
}

Eigen::Vector3d learnSPH::Kernel::kernelGrad(Eigen::Vector3d x_i,
                                             Eigen::Vector3d x_j, double h)
{
    Eigen::Vector3d dist = (x_i - x_j);
    Eigen::Vector3d diff = 2*dist;
    for (int i = 0; i < 3; i++)
        diff[i] = 1/sqrt(diff[i]);

    double q = dist.norm() / h;
    return pow(h, -3) * cubicSplineGrad(q) * 0.5*h * diff;
}

double learnSPH::Kernel::cubicSpline(const double q)
{
    assert(q >= 0.0);
        
    constexpr double alpha = 3.0 / (2.0 * PI);
    double value = 0.0;
        
    if (q < 1.0) {
        value = (2.0/3.0 - q*q + 0.5*q*q*q);
    }
    else if (q < 2.0) {
        value = pow((2.0-q), 3) / 6.0;
    }
        
    return alpha * value;
}

double learnSPH::Kernel::cubicSplineGrad(const double q)
{
    assert(q >= 0.0);
        
    constexpr double alpha = 3.0 / (2.0 * PI);
    double value = 0.0;
        
    if (q < 1.0) {
        value = (3.0/2.0)*q*q - 2*q;
    }
    else if (q < 2.0) {
        value = -(1.0/2.0)*(2-q)*(2-q);
    }
        
    return alpha * value;
}

