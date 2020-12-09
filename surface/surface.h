#include <Eigen/Dense>
#include <functional>
#include <vector>

namespace learnSPH::Surface {
    //void printCudaVersion();
    
    std::vector<Eigen::Vector3d> marchCubes(
        Eigen::Vector3i volDim,
        std::vector<Eigen::Vector3f> &volVerts,
        std::vector<float> &volSDF);
    

    void discretizeSDF(
        const Eigen::Vector3f volSize,
        const Eigen::Vector3i volDims,
        std::function<float(Eigen::Vector3f)> const &sdf,
        std::vector<Eigen::Vector3f> *pOutVolVerts,
        std::vector<float> *pOutVolSDF);
}
