#include <Eigen/Dense>
#include <functional>
#include <vector>

typedef Eigen::Matrix<size_t, 3, 1> IndexVector3d;
typedef Eigen::Matrix<size_t, 3, 1> Dimensions3d;  

namespace learnSPH::Surface {
    //void printCudaVersion();
    
    std::vector<Eigen::Vector3d> marchCubes(
        std::vector<Eigen::Vector3d> const &gridVerts,
        std::vector<double> const &gridSDF,
        const Dimensions3d gridDim);
    

    void discretizeSDF(
        const Dimensions3d gridDim,
        const Eigen::Vector3d gridSize,
        std::function<double(Eigen::Vector3d)> const &sdf,
        std::vector<Eigen::Vector3d> *pOutGridVerts,
        std::vector<double> *pOutGridSDF);
}
