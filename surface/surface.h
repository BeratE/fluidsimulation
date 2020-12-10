#include <Eigen/Dense>
#include <functional>
#include <vector>

namespace learnSPH::Surface {
    //void printCudaVersion();
    
    void marchCubes(
        const Eigen::Vector3i volDim,
        const std::vector<float> &volSDF,
        const std::vector<Eigen::Vector3f> &volVerts,
        std::vector<Eigen::Vector3f> &outVertices,
        std::vector<std::array<int, 3>> &outTriangles);
    

    void discretizeSDF(
        const Eigen::Vector3f volSize,
        const Eigen::Vector3i volDims,
        std::function<float(Eigen::Vector3f)> const &sdf,
        std::vector<float> *pOutVolSDF,
        std::vector<Eigen::Vector3f> *pOutVolVerts);
}
