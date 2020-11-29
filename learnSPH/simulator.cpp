#include <Eigen/Dense>
#include "simulator.h"
#include "kernel.h"
#include <math.h>
#include "util/vtk_writer.h"
#include "config.h"


// using namespace learnSPH;
// using namespace learnSPH::Simulator;
// using namespace learnSPH::ParticleSystem;

// void Simulator::semiImplicitEuler(ParticleSystem::FluidSystem& fluid,
//                                   const CompactNSearch::NeighborhoodSearch& nsearch,
//                                   const double defaultTimeStep,
//                                   const double epsilon,
//                                   bool smoothing)
// {
//     // smoothin length 
//     const double h = Kernel::Parameter::TUNING * fluid.particleRadius * 2;
//     const double timeStepCFL = fluid.getTimeCFL();
//     const double timeStep = std::min(defaultTimeStep, timeStepCFL);

//     // get neighborhood information of fluid particle point set
//     CompactNSearch::PointSet const& fluidPS = nsearch.point_set(fluid.id);
//     fluid.positions.resize(fluidPS.n_points());
//     fluid.velocities.resize(fluidPS.n_points());
//     fluid.accelerations.resize(fluidPS.n_points());

//     for (size_t fpI = 0; fpI < fluidPS.n_points(); fpI++) {
//         fluid.velocities[fpI] += timeStep * fluid.accelerations[fpI];
//     }

//     // iterate fluid particles 
//     for (size_t fpI = 0; fpI < fluidPS.n_points(); fpI++) {
//         const Eigen::Vector3d &fpPos = fluid.positions[fpI];
//         const Eigen::Vector3d &fpVel = fluid.velocities[fpI];
//         Eigen::Vector3d fpVelStar = fpVel;

//         if (smoothing) {
//             Eigen::Vector3d fpVelSumOverNeighbors = Eigen::Vector3d::Zero();
//             /* No need to pay attention to the particle itself when calculating fpSmoothedVelSum,
//              * because vj - vi == 0 if i == j */
//             for (size_t fpN = 0; fpN < fluidPS.n_neighbors(fluid.id, fpI); fpN++) {
//                 const unsigned int fnI = fluidPS.neighbor(fluid.id, fpI, fpN);
//                 fpVelSumOverNeighbors +=
//                     (fluid.velocities[fnI] - fluid.velocities[fpI]) /
//                     (fluid.densities[fpI] + fluid.densities[fnI]) *
//                     Kernel::CubicSpline::weight(fpPos, fluid.positions[fnI], h);
//             }

//             fpVelSumOverNeighbors *= 2 * fluid.particleMass * epsilon;
//             fpVelStar += fpVelSumOverNeighbors;
//         }
//         fluid.positions[fpI] += timeStep * fpVelStar;
//     }
// }

// void Simulator::simulate(ParticleSystem::FluidSystem& fluid,
//                          std::vector<ParticleSystem::BoundarySystem>& boundaries,
//                          const double defaultTimeStep,
//                          const int simulationSteps,
//                          const int zSortSkip,
//                          const std::string fileBaseName)
// {
//     for (const BoundarySystem& boundary : boundaries) {
//         int boundaryIdx = 0;
//         std::vector<Eigen::Vector3d> boundaryPositions = boundary.positions;
//         save_particles_to_vtk(SOURCE_DIR + std::string("/res/simulation/")
//                               + fileBaseName + std::string("_boundary")
//                               + std::to_string(boundaryIdx)
//                               + std::string(".vtk"), boundaryPositions);
//         boundaryIdx++;
//     }
	
//     std::vector<Eigen::Vector3d> prevPositions = fluid.positions;
//     save_particles_to_vtk(SOURCE_DIR + std::string("/res/simulation/") + fileBaseName + std::string("0.vtk"), prevPositions);
//     double timeStep = std::min(defaultTimeStep, fluid.getTimeCFL());
//     double nextStopTime = defaultTimeStep;
//     double simulatedTime = 0.0;
//     double prevTime = 0.0;
//     int fileNr = 1;
//     int simulationStep = 0;

//     const double h = 2.0 * fluid.particleRadius * Kernel::Parameter::TUNING;
//     CompactNSearch::NeighborhoodSearch nsearch(Kernel::CubicSpline::support(h));
    
//     fluid.id = nsearch.add_point_set(fluid.positions.front().data(), fluid.positions.size());
    
//     for (BoundarySystem& boundary : boundaries) {
//         boundary.id = nsearch.add_point_set(boundary.positions.front().data(), boundary.positions.size());
//     }

//     nsearch.find_neighbors();
//     estimateFluidDensity(fluid, nsearch);

//     while (fileNr <= simulationSteps) {
//         std::cout << simulatedTime << std::endl;
//         if (simulationStep % zSortSkip == 0) {
//             nsearch.z_sort();
//             auto const& fluidPS = nsearch.point_set(fluid.id);
//             fluidPS.sort_field(fluid.positions.data());
//             fluidPS.sort_field(fluid.densities.data());
//             fluidPS.sort_field(fluid.velocities.data());
//             fluidPS.sort_field(fluid.accelerations.data());
//             for (BoundarySystem& boundary : boundaries) {
//                 auto const& boundaryPS = nsearch.point_set(boundary.id);
//                 boundaryPS.sort_field(boundary.positions.data());
//                 boundaryPS.sort_field(boundary.volumes.data());
//             }
//         }
//         nsearch.find_neighbors();
//         estimateFluidDensity(fluid, nsearch);
//         fluid.updateAccelerations(boundaries, nsearch);
//         semiImplicitEuler(fluid, nsearch, defaultTimeStep);
//         simulatedTime += timeStep;
//         if (simulatedTime >= nextStopTime) {
//             std::string filename = SOURCE_DIR + std::string("/res/simulation/") + fileBaseName + std::to_string(fileNr) + std::string(".vtk");
//             std::vector<Eigen::Vector3d> interpolatedPositions = interpolatePositions(prevPositions, fluid.positions, prevTime, simulatedTime, nextStopTime);
//             save_particles_to_vtk(filename, interpolatedPositions);
//             fileNr++;
//             nextStopTime += defaultTimeStep;
//         }
//         prevPositions = fluid.positions;
//         prevTime += timeStep;

//         timeStep = std::min(defaultTimeStep, fluid.getTimeCFL());
//         simulationStep++;
//     }
// }


