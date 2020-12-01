#include "catch.hpp"
#include <Eigen/Dense>
#include <tuple>
#include <stdio.h>
#include "kernel.h"
#include <string>
#include "emitter.h"
#include "particlesystem.h"
#include "fluidsystem.h"
#include "util/vtk_writer.h"
#include "util/config.h"

using namespace learnSPH;
using namespace learnSPH::Kernel;

TEST_CASE("FluidBox", "Generate fluid particles in a box")
{
    std::cout << "Generating Fluid Particles in a Box.." << std::endl;

    const double particleDiameter = 0.1;
    FluidSystem particles = Emitter().sampleFluidBox(Eigen::Vector3d(0, 0, 0),
                                                     Eigen::Vector3d(1, 1, 1),
                                                     particleDiameter);
    
    const double r = CubicSpline::support(particles.smoothingLength());
    CompactNSearch::NeighborhoodSearch nsearch(r);
    particles.addToNeighborhood(nsearch);
    
    nsearch.find_neighbors();
    
    particles.estimateDensity(nsearch);

    std::stringstream filename;
    filename << SOURCE_DIR << "/res/test_fluidbox.vtk";
    save_particles_to_vtk(filename.str(), particles.getPositions(), particles.getDensities());

    std::cout << "Results saved to " << filename.str() << std::endl;
}


TEST_CASE("FluidBoxBoundary", "Generate fluid particles enclosed in a boundary box")
{
    std::cout << "Generating Fluid Particles in a Boundary (hollow) box.." << std::endl;
    
    const double particleDiameter = 0.1;

    FluidSystem particles = Emitter().sampleFluidBox(Eigen::Vector3d(0, 0, 0),
                                                     Eigen::Vector3d(1, 1, 0.6),
                                                     particleDiameter);

    BoundarySystem boundary = Emitter().sampleBoundaryHollowBox(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                                Eigen::Vector3d(1.0, 1.0, 1.0),
                                                                particleDiameter);
    std::vector<BoundarySystem> boundaries;
    
    boundaries.push_back(boundary);
    
    // Compute neighborhood information of fluid particles
    const double r = Kernel::CubicSpline::support(particles.smoothingLength());
    CompactNSearch::NeighborhoodSearch nsearch(r);
    particles.addToNeighborhood(nsearch);
    boundaries[0].addToNeighborhood(nsearch);
    
    nsearch.find_neighbors();
    
    particles.estimateDensity(nsearch, boundaries);

    std::stringstream f1;
    f1 << SOURCE_DIR << "/res/test_fluid_in_box.vtk";
    save_particles_to_vtk(f1.str(), particles.getPositions(), particles.getDensities());
    std::cout << "Results saved to " << f1.str() << std::endl;
    
    std::stringstream f2;
    f2 << SOURCE_DIR  << "/res/test_boundary_box.vtk";
    save_particles_to_vtk(f2.str(), boundaries[0].getPositions(), boundaries[0].getVolumes());
    std::cout << "Results saved to " << f2.str() << std::endl;
}
