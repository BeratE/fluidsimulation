#include "catch.hpp"
#include <Eigen/Dense>
#include <tuple>
#include <stdio.h>
#include "kernel.h"
#include <string>
#include "particlesystem.h"
#include "util/vtk_writer.h"
#include "config.h"

using namespace learnSPH;
using namespace learnSPH::ParticleSystem;

TEST_CASE("Generate fluid particles in a box", "[box]")
{
    std::cout << "Generating box sample scene.." << std::endl;

    // Generate particles
    const double restDensity = 1.0;
    const double particleDiameter = 0.1;
    const double smoothingLength = particleDiameter * Kernel::Parameter::TUNING;
    // Sample Particles in a Box
    FluidSystem particles
        = sampleFluidBox(Eigen::Vector3d(0, 0, 0),
                         Eigen::Vector3d(1, 1, 1),
                         restDensity, particleDiameter);
    // Compute neighborhood information of fluid particles
    CompactNSearch::NeighborhoodSearch nsearch(Kernel::CubicSpline::support(smoothingLength));
    particles.id = nsearch.add_point_set(particles.positions.front().data(),
                                         particles.positions.size());
    nsearch.find_neighbors();
    
    estimateFluidDensity(particles, nsearch);

    std::stringstream filename;
    filename << SOURCE_DIR << "/res/test_p1_d01_n12.vtk";
    save_particles_to_vtk(filename.str(), particles.positions, particles.densities);

    std::cout << "completed!" << std::endl;
    std::cout << "The scene files have been saved in the folder "
        "`<source_folder>/res`. You can visualize them with Paraview."
              << std::endl;
}

TEST_CASE("Sample positions in a triangle", "[pos_triangle]")
{
    double samplingDistance = 0.1;
    std::vector<Eigen::Vector3d> positions = samplePositionsTriangle(Eigen::Vector3d(0.0,  0.0,  0.0),
                                                                     Eigen::Vector3d(0.06, 0.0,  0.0),
                                                                     Eigen::Vector3d(0.06, 0.06, 0.0),
                                                                     samplingDistance);
    std::stringstream filename;
    filename << SOURCE_DIR << "/res/test_tiangle.vtk";
    save_particles_to_vtk(filename.str(), positions);

    std::cout << "completed!" << std::endl;
    std::cout << "The scene files have been saved in the folder "
        "`<source_folder>/res`. You can visualize them with Paraview."
              << std::endl;
}

TEST_CASE("Sample hollow boundary box", "[pos_box]")
{
    double samplingDistance = 0.1;
    std::vector<Eigen::Vector3d> positions = samplePositionsBox(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                                Eigen::Vector3d(1.0, 1.0, 1.0),
                                                                samplingDistance);
    std::stringstream filename;
    filename << SOURCE_DIR << "/res/test_box.vtk";
    save_particles_to_vtk(filename.str(), positions);

    std::cout << "completed!" << std::endl;
    std::cout << "The scene files have been saved in the folder "
        "`<source_folder>/res`. You can visualize them with Paraview."
              << std::endl;
}
