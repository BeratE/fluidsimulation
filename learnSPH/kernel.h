#pragma once
#include <cassert>
#include <Eigen/Dense>
#include <CompactNSearch/CompactNSearch.h>

namespace learnSPH::Kernel {
    namespace Parameter {
        constexpr double TUNING = 1.4;
        constexpr double EPSILON = 0.5;
    } // namespace Parameter

    // Cubic Spline Kernel
    namespace CubicSpline {
        double support(double h);
        double weight(Eigen::Vector3d x_i, Eigen::Vector3d x_j, double h);
        Eigen::Vector3d gradWeight(Eigen::Vector3d x_i, Eigen::Vector3d x_j, double h);
        double cubicSpline(const double q);
        double gradCubicSpline(const double q);

        class Table {
        public:
            Table() {}
            Table(double smoothingLength, size_t numBins = 100);
            void generateTable(double smoothingLength, size_t numBins = 100);

            double weight(Eigen::Vector3d x_i, Eigen::Vector3d x_j);
            Eigen::Vector3d gradWeight(Eigen::Vector3d x_i, Eigen::Vector3d x_j);

        private:
            bool m_isInit = false;
            double m_support;
            double m_stepSize;
            std::vector<double> m_weights;
            std::vector<Eigen::Vector3d> m_weightGrads;
        };
    } // namespace CubicSpline
    
} // namespace learnSPH::Kernel
