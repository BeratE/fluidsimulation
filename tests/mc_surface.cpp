#include "catch.hpp"
#include <tuple>
#include <string>
#include <iostream>
#include <Eigen/Dense>
#include "vtk_writer.h"
#include "config.h"
#include "surface/surface.h"
#include <omp.h>

using namespace learnSPH;
using namespace learnSPH::Surface;

TEST_CASE("Construction", "") {
    //printCudaVersion();
    
    const double start_t = omp_get_wtime();

    const Dimensions3d gridDim(30, 30, 30);
    const Eigen::Vector3d gridSize(1.0, 1.0, 1.0);
    
    std::vector<Eigen::Vector3d> gridVerts;
    std::vector<double> gridSDF;
    
    SECTION("SphereSDF") {
        const double radius = 0.2;
        const Eigen::Vector3d origin(0.5, 0.5, 0.5);
        auto sdf = [origin,radius](Eigen::Vector3d x)
            { return (x-origin).norm() - radius; };
        
        discretizeSDF(gridDim, gridSize, sdf,
                  &gridVerts, &gridSDF);
    }
    SECTION("TorusSDF") {
        const double r = 0.1;
        const double R = 0.3;
        const Eigen::Vector3d origin(0.5, 0.5, 0.5);
        auto sdf = [r, R, origin](Eigen::Vector3d x)
            {
                x = x - origin;
                return r*r-pow(sqrt(x(0)*x(0)+x(1)*x(1))-R, 2)-x(2)*x(2);
            };
        
        discretizeSDF(gridDim, gridSize, sdf,
                  &gridVerts, &gridSDF);
    }

    std::vector<Eigen::Vector3d> triangles =
        marchCubes(gridVerts, gridSDF, gridDim);

    std::stringstream filename;
    filename << SOURCE_DIR << "/res/surface/simple_surface.vtk";
        
    save_particles_to_vtk(filename.str(), triangles);

    const double end_t = omp_get_wtime();
    const double delta_t = end_t - start_t;

    std::cout << "Runtime: " << delta_t << std::endl;
}
