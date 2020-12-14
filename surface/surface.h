#include <Eigen/Dense>
#include <functional>
#include <vector>
#include <learnSPH/system/fluidsystem.h>

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
        const learnSPH::System::FluidSystem &system, 
        const double c, 
        const double samplingDistance,
        std::vector<double> *pOutVolSDF,
        std::vector<Eigen::Vector3d> *pOutVolVerts,
        Eigen::Vector3i* pOutDims);

    size_t getVertIdx(Eigen::Vector3i pos, Eigen::Vector3i volDim);

}
