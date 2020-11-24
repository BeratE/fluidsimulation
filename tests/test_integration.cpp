#include "catch.hpp"
#include "particlesystem.h"
#include "simulator.h"
#include "kernel.h"
#include "util/vtk_writer.h"
#include "config.h"

using namespace learnSPH;
using namespace learnSPH::ParticleSystem;
using namespace learnSPH::Simulator;

TEST_CASE("Test semi-implicit Euler for a cube of particles with gravity only", "[semiEuler]") {
	std::cout << "Testing semi-implicit Euler.." << std::endl;

    // Generate particles
    const double restDensity = 1.0;
    const double particleDiameter = 0.1;
    const double smoothingLength = particleDiameter * Kernel::Parameter::TUNING;

    // Sample Particles in a Box
    FluidSystem particles
        = sampleFluidBox(Eigen::Vector3d(0, 0, 0),
            Eigen::Vector3d(1, 1, 1),
            restDensity, particleDiameter);

    particles.setAccelerationsToGravity();

    // Compute neighborhood information of fluid particles
    CompactNSearch::NeighborhoodSearch nsearch(Kernel::CubicSpline::support(smoothingLength));
    particles.id = nsearch.add_point_set(particles.positions.front().data(),
        particles.positions.size());

    std::string filename = SOURCE_DIR + std::string("/res/integration/test_semi_euler") + std::to_string(0) + std::string(".vtk");
    
    save_particles_to_vtk(filename, particles.positions);
    for (size_t steps = 1; steps < 10; steps++) {
        filename = SOURCE_DIR + std::string("/res/integration/test_semi_euler") + std::to_string(steps) + std::string(".vtk");
        nsearch.find_neighbors();
        estimateFluidDensity(particles, nsearch);
        semiImplicitEuler(particles, 0.1, nsearch, 0.5);
        save_particles_to_vtk(filename, particles.positions);
    }

    std::cout << "completed" << std::endl;
    std::cout << "The scene files have been saved in the folder `<source_folder>/res/integration`. You can visualize them with Paraview." << std::endl;
}
