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

    const Eigen::Vector3i numVerts(30, 30, 30);
    const Eigen::Vector3f gridSize(1.0, 1.0, 1.0);
    
    std::vector<Eigen::Vector3f> gridVerts;
    std::vector<float> gridSDF;
    
    SECTION("SphereSDF") {
        const float radius = 0.2;
        const Eigen::Vector3f origin(0.5, 0.5, 0.5);
        auto sdf = [origin,radius](Eigen::Vector3f x)
            { return (x-origin).norm() - radius; };
        
        discretizeSDF(gridSize, numVerts, sdf,
                      &gridSDF, &gridVerts);
    }
    SECTION("TorusSDF") {
        const float r = 0.1;
        const float R = 0.3;
        const Eigen::Vector3f origin(0.5, 0.5, 0.5);
        auto sdf = [r, R, origin](Eigen::Vector3f x)
            {
                x = x - origin;
                return r*r-pow(sqrt(x(0)*x(0)+x(1)*x(1))-R, 2)-x(2)*x(2);
            };
        
        discretizeSDF(gridSize, numVerts, sdf,
                      &gridSDF, &gridVerts);
    }

    std::vector<Eigen::Vector3f> vertices;
    std::vector<std::array<int, 3>> triangles;
    marchCubes(numVerts, gridSDF, gridVerts,
               vertices, triangles);

    std::stringstream filename;
    filename << SOURCE_DIR << "/res/surface/simple_surface.vtk";

    std::vector<Eigen::Vector3d> verts;
    for (auto &v : vertices) {
        verts.push_back(v.cast<double>());
    }

    save_mesh_to_vtk(filename.str(), vertices, triangles);
    //save_particles_to_vtk(filename.str(), verts);

    const double end_t = omp_get_wtime();
    const double delta_t = end_t - start_t;

    std::cout << "Runtime: " << delta_t << std::endl;
}
