#pragma once
#include <Eigen/Dense>
#include <functional>
#include <vector>
#include <learnSPH/system/fluidsystem.h>

namespace learnSPH::Surface {
    //void printCudaVersion();
    
    class SurfaceInformation {
    public:
        SurfaceInformation(
            const std::vector<Eigen::Vector3d> positions,
            const std::vector<double> normalizedDensities,
            const Kernel::CubicSpline::Table kernelLookup,
            double smoothingLength,
            std::string filename);

        // Setter & Getter
        const std::vector<Eigen::Vector3d> getPositions() { return m_positions; }
        const std::vector<double> getNormalizedDensities() { return m_normalizedDensities; }
        const Kernel::CubicSpline::Table getKernelLookup() { return m_kernelLookup; };
        double getSmoothingLength() { return m_smoothingLength; };
        std::string getFilename() { return m_filename;  };

    protected:
        const std::vector<Eigen::Vector3d> m_positions;
        const std::vector<double> m_normalizedDensities;
        const Kernel::CubicSpline::Table m_kernelLookup;
        double m_smoothingLength;
        std::string m_filename;
    };

    void marchCubes(
        const Eigen::Vector3i volDim,
        const std::vector<double> &volSDF,
        const std::vector<Eigen::Vector3d> &volVerts,
        std::vector<Eigen::Vector3d> &outVertices,
        std::vector<std::array<int, 3>> &outTriangles);
    

    void discretizeSDF(
        const Eigen::Vector3d volSize,
        const Eigen::Vector3i volDims,
        std::function<double(Eigen::Vector3d)> const &sdf,
        std::vector<double> *pOutVolSDF,
        std::vector<Eigen::Vector3d> *pOutVolVerts);

    void discretizeFluidSystemSDF(
        const std::vector<Eigen::Vector3d>& positions,
        const std::vector<double>& normalizedDensities,
        const Kernel::CubicSpline::Table& kernelLookup,
        double smoothingLength,
        const double c,
        const double samplingDistance,
        std::vector<double>* pOutVolSDF,
        std::vector<Eigen::Vector3d>* pOutVolVerts,
        Eigen::Vector3i* pOutDims);

    size_t getVertIdx(Eigen::Vector3i pos, Eigen::Vector3i volDim);

    void boundingBox(const std::vector<Eigen::Vector3d>& positions, Eigen::Vector3d& bottomLeft, Eigen::Vector3d& upperRight);
}
