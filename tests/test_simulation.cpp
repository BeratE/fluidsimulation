#include "catch.hpp"
#include "particlesystem.h"
#include "simulator.h"
#include "kernel.h"
#include "util/vtk_writer.h"
#include "config.h"

using namespace learnSPH;

// TEST_CASE("semiEuler", "Test semi-implicit Euler for a cube of particles with gravity only") {
//     std::cout << "Testing semi-implicit Euler.." << std::endl;

//     // Generate particles
//     const double restDensity = 1.0;
//     const double particleDiameter = 0.1;
//     const double smoothingLength = particleDiameter * Kernel::Parameter::TUNING;

//     // Sample Particles in a Box
//     FluidSystem particles = sampleFluidBox(Eigen::Vector3d(0, 0, 0),
//                                            Eigen::Vector3d(1, 1, 1),
//                                            restDensity, particleDiameter);

//     particles.setAccelerationsToGravity();

//     // Compute neighborhood information of fluid particles
//     CompactNSearch::NeighborhoodSearch nsearch(Kernel::CubicSpline::support(smoothingLength));
//     particles.id = nsearch.add_point_set(particles.positions.front().data(),
//                                          particles.positions.size());

//     std::string filename = SOURCE_DIR + std::string("/res/integration/test_semi_euler_smooth_on")
//         + std::to_string(0) + std::string(".vtk");
//     save_particles_to_vtk(filename, particles.positions);

//     for (size_t steps = 1; steps < 100; steps++) {
//         filename = SOURCE_DIR +
//             std::string("/res/integration/test_semi_euler_smooth_on") +
//             std::to_string(steps) + std::string(".vtk");
//         nsearch.find_neighbors();
//         estimateFluidDensity(particles, nsearch);
//         semiImplicitEuler(particles, nsearch);
//         save_particles_to_vtk(filename, particles.positions);
//     }

//     filename = SOURCE_DIR + std::string("/res/integration/test_semi_euler_smooth_off")
//         + std::to_string(0) + std::string(".vtk");
//     save_particles_to_vtk(filename, particles.positions);

//     for (size_t steps = 1; steps < 100; steps++) {
//         filename = SOURCE_DIR +
//             std::string("/res/integration/test_semi_euler_smooth_off") +
//             std::to_string(steps) + std::string(".vtk");
//         nsearch.find_neighbors();
//         estimateFluidDensity(particles, nsearch);
//         semiImplicitEuler(particles, nsearch, 0.01, Kernel::Parameter::EPSILON, false);
//         save_particles_to_vtk(filename, particles.positions);
//     }

//     std::cout << "completed" << std::endl;
//     std::cout << "The scene files have been saved in the folder `<source_folder>/res/integration`. You can visualize them with Paraview." << std::endl;
// }

// TEST_CASE("Test time constant simulation for a cube of particles with gravity in a box", "[simulation]") {
//     std::cout << "Testing time constant simulation.." << std::endl;

//     // Generate particles
//     const double restDensity = 1000;
//     const double particleDiameter = 0.2;
//     const double smoothingLength = particleDiameter * Kernel::Parameter::TUNING;

//     // Sample Particles in a Box
//     FluidSystem particles
//         = sampleFluidBox(Eigen::Vector3d(0, 0, 0),
//             Eigen::Vector3d(1, 1, 1),
//             restDensity, particleDiameter);

//     particles.setAccelerationsToGravity();
//     particles.viscosity = 0.1;
//     std::vector<BoundarySystem> boundaries;
//     std::vector<Eigen::Vector3d> boundaryPositions = samplePositionsBox(Eigen::Vector3d(-0.1, -0.1, -0.1), Eigen::Vector3d(2, 2, 2), particleDiameter * 0.5);
//     BoundarySystem boundary = createBoundary(boundaryPositions, 2600, 1000, particleDiameter);
//     boundary.viscosity = 1;
//     boundaries.push_back(boundary);
    
//     simulate(particles, boundaries, 0.02, 10, 50, std::string("time_consistent_simulation"));
// }
