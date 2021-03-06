#pragma once
#include <Eigen/Dense>
#include <functional>
#include <vector>
#include <learnSPH/system/fluidsystem.h>
#include <common/inc/vtk_writer.h>

namespace learnSPH::Surface {
    //void printCudaVersion();
    
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
        double samplingDistance,
        double paramC,
        std::vector<double>* pOutVolSDF,
        std::vector<Eigen::Vector3d>* pOutVolVerts,
        Eigen::Vector3i* pOutDims);   
    
    size_t getVertIdx(Eigen::Vector3i pos, Eigen::Vector3i volDim);

    void boundingBox(const std::vector<Eigen::Vector3d>& positions,
                     Eigen::Vector3d& bottomLeft, Eigen::Vector3d& upperRight);
}
