#include "kernel.h"

double learnSPH::Kernel::kernel(Eigen::Vector3d x_i, Eigen::Vector3d x_j,
                                double h)
{
    double q = (x_j - x_i).norm() / h;

    return (1/pow(h, 3)) * cubicSpline(q);}

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
	// ...
	return 0.0;
}

