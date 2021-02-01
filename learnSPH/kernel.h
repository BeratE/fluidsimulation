#pragma once

#include <cassert>
#include <Eigen/Dense>
#include <CompactNSearch/CompactNSearch.h>
#define _USE_MATH_DEFINES
#include <math.h>

namespace learnSPH::Kernel {
    namespace Parameter {
        constexpr double TUNING = 1.4;
        constexpr double EPSILON = 0.5;
    } // namespace Parameter

    // Cubic Spline Kernel
    namespace CubicSpline {
        double support(double h);
        double weight(Eigen::Vector3d x_i,
                      Eigen::Vector3d x_j, double h);
        Eigen::Vector3d gradWeight(Eigen::Vector3d x_i,
                                   Eigen::Vector3d x_j, double h);
        
        double cubicSpline(const double q);
        double gradCubicSpline(const double q);

        class Table {
        public:
            Table() {};
            Table(double smoothingLength, size_t numBins = 1000);
            void generateTable(double smoothingLength,
                               size_t numBins = 1000);

            double weight(Eigen::Vector3d x_i,
                          Eigen::Vector3d x_j) const;
            Eigen::Vector3d gradWeight(Eigen::Vector3d x_i,
                                       Eigen::Vector3d x_j) const;

        private:
            bool m_isInit = false;
            double m_support;
            double m_stepSize;
            double m_smoothingLength;
            std::vector<double> m_weights;
            std::vector<double> m_gradMagnitudes;
        };
    } // namespace CubicSpline

    namespace Cohesion {
        double weight(const double r, const double c);
        class Table {
        public:
            Table() {};
            Table(const double c, size_t numBins = 1000);
            void generateTable(const double c,
                size_t numBins = 1000);
            double weight(const double r);
        private:
            std::vector<double> m_weights;
            double m_support = 2.0;
            double m_stepSize;
        };
    }

    namespace Adhesion {
        double weight(const double r, const double c);
        class Table {
        public:
            Table() {};
            Table(const double c, size_t numBins = 1000);
            void generateTable(const double c,
                size_t numBins = 1000);
            double weight(const double r);
        private:
            std::vector<double> m_weights;
            double m_support;
            double m_stepSize;
            double m_c;
        };
    }
} // namespace learnSPH::Kernel
