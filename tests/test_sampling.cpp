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

TEST_CASE("FluidBox", "Generate fluid particles in a box")
{
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
}


TEST_CASE("FluidBoxBoundary", "Generate fluid particles enclosed in a boundary box")
{
    // Generate particles
    const double restDensity = 1.0;
    const double particleDiameter = 0.1;
    const double smoothingLength = particleDiameter * Kernel::Parameter::TUNING;
    // Sample Particles in a Box
    FluidSystem particles = sampleFluidBox(Eigen::Vector3d(0, 0, 0),
                                           Eigen::Vector3d(1, 1, 0.6),
                                           restDensity, particleDiameter);
    
    std::vector<Eigen::Vector3d> boxpos = samplePositionsBox(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                             Eigen::Vector3d(1.0, 1.0, 1.0),
                                                             particleDiameter * 0.9);
    std::vector<BoundarySystem> boundaries;
    
    boundaries.push_back(createBoundary(boxpos, 2.6, 1, particleDiameter));
    
    // Compute neighborhood information of fluid particles
    CompactNSearch::NeighborhoodSearch nsearch(Kernel::CubicSpline::support(smoothingLength));
    particles.id = nsearch.add_point_set(particles.positions.front().data(),
                                         particles.positions.size());
    boundaries[0].id = nsearch.add_point_set(boundaries[0].positions.front().data(),
                                        boundaries[0].positions.size());
    
    nsearch.find_neighbors();
    
    estimateFluidDensity(particles, boundaries, nsearch);

    std::stringstream f1;
    f1 << SOURCE_DIR << "/res/test_fluid_in_box.vtk";
    save_particles_to_vtk(f1.str(), particles.positions, particles.densities);
    std::stringstream f2;
    f2 << SOURCE_DIR  << "/res/test_boundary_box.vtk";
    save_particles_to_vtk(f2.str(), boundaries[0].positions);

}

TEST_CASE("SampleTriangle", "Sample positions in a triangle")
{
    double samplingDistance = 0.1;
    std::vector<Eigen::Vector3d> positions = samplePositionsTriangle(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                                     Eigen::Vector3d(0.5, 0.0, 0.0),
                                                                     Eigen::Vector3d(1.0, 1.0, 1.0),
                                                                     samplingDistance);
    std::stringstream filename;
    filename << SOURCE_DIR << "/res/test_tiangle.vtk";
    save_particles_to_vtk(filename.str(), positions);
}

TEST_CASE("SampleBox", "Sample positions on the outer shell of a box")
{
    double samplingDistance = 0.1;
    std::vector<Eigen::Vector3d> positions = samplePositionsBox(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                                Eigen::Vector3d(1.0, 1.0, 1.0),
                                                                samplingDistance);
    std::stringstream filename;
    filename << SOURCE_DIR << "/res/test_box.vtk";
    save_particles_to_vtk(filename.str(), positions);
}
