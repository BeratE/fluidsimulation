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




double Cohesion::weight(const double r, const double c) {
    const double alpha = 32.0 / (M_PI * c * c * c * c * c * c * c * c * c);
    if (0.0 <= r && r <= c / 2.0) {
        return (alpha * 2 * (c - r) * (c - r) * (c - r) * r * r * r) - ((c * c * c * c * c * c) / 64.0);
    }
    else if (c / 2.0 < r && r <= c) {
        return alpha * (c - r) * (c - r) * (c - r) * r * r * r;
    }
    return 0.0;
}

double Adhesion::weight(const double r, const double c) {
    const double alpha = 0.007 / pow(c, 3.25);

    if (c / 2.0 <= r && r <= c) {
        return alpha * pow((-4.0 * r * r) / c + 6.0 * r - 2.0 * c, (1.0 / 4.0));
    }
    return 0;
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


Cohesion::Table::Table(const double c, size_t numBins)
{
    generateTable(c, numBins);
}

void Cohesion::Table::generateTable(double c, size_t numBins)
{
    m_stepSize = m_support / numBins;

    for (int i = 0; i < numBins; i++) {
        double r = (float(i) / float(numBins)) * m_support;
        m_weights.push_back(Cohesion::weight(r, c));
    }
}

double Cohesion::Table::weight(const double r) {
    if (r > m_support) {
        return 0.0;
    }

    size_t i = (size_t)floor(r / m_stepSize);
    return m_weights[i];
}

Adhesion::Table::Table(const double c, size_t numBins)
{
    generateTable(c, numBins);
}

void Adhesion::Table::generateTable(double c, size_t numBins)
{
    m_c = c;
    m_support = c / 2.0;
    m_stepSize = m_support / numBins;

    for (int i = 0; i < numBins; i++) {
        double r = c / 2.0 + i * m_stepSize;
        m_weights.push_back(Adhesion::weight(r, c));
    }
}

double Adhesion::Table::weight(const double r) {
    if (r - m_c / 2.0 > m_support) {
        return 0.0;
    }

    size_t i = (size_t)floor( ( r - m_c / 2.0 ) / m_stepSize);
    return m_weights[i];
}
