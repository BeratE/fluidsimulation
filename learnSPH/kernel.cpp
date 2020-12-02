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
    double q = (x_i - x_j).norm() / h;
    return pow(h, -3) * cubicSpline(q);
}

Eigen::Vector3d CubicSpline::gradWeight(Eigen::Vector3d x_i, Eigen::Vector3d x_j, double h)
{
    Eigen::Vector3d dist = (x_i - x_j);
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


// Table
CubicSpline::Table::Table(double smoothingLength, size_t numBins)
{
    generateTable(smoothingLength, numBins);
}

void CubicSpline::Table::generateTable(double smoothingLength, size_t numBins)
{
    m_smoothingLength = smoothingLength;
    m_support = CubicSpline::support(smoothingLength);
    m_stepSize = m_support/numBins;

    for (int i = 0; i < numBins; i++) {
        double d = i*m_stepSize;
        double q = d /smoothingLength;
        m_weights.push_back(CubicSpline::weight(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                Eigen::Vector3d(d, 0.0, 0.0),
                                                smoothingLength));
        m_gradMagnitudes.push_back(CubicSpline::gradCubicSpline(q)
                                   * pow(m_smoothingLength, -4));
    }
    m_isInit = true;
}

double CubicSpline::Table::weight(Eigen::Vector3d x_i,
                                  Eigen::Vector3d x_j) const
{
    assert(m_isInit);
    
    double d = (x_i - x_j).norm();
    if (d > m_support)
        return 0.0;
    
    size_t i = (size_t)floor(d/m_stepSize);
    return m_weights[i];
}

Eigen::Vector3d CubicSpline::Table::gradWeight(Eigen::Vector3d x_i,
                                               Eigen::Vector3d x_j) const
{
    assert(m_isInit);
    
    Eigen::Vector3d posDiff = (x_i - x_j);
    double d = posDiff.norm();
    if (d > m_support)
        return Eigen::Vector3d(0.0, 0.0, 0.0);
    
    size_t i = (size_t)floor(d/m_stepSize);
    return m_gradMagnitudes[i] * posDiff.normalized();
}
