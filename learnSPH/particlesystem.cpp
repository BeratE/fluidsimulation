#include "particlesystem.h"
#include "kernel.h"

using namespace learnSPH;
using namespace CompactNSearch;

double ParticleSystem::smoothingLength() const
{
    return 2 * m_particleRadius * Kernel::Parameter::TUNING;
}

double ParticleSystem::particleMass() const
{
    return pow(2 * m_particleRadius, 3) * m_restDensity;
}

void ParticleSystem::addToNeighborhood(NeighborhoodSearch &nsearch)
{
    m_pointSetID = nsearch.add_point_set(m_positions.front().data(),
                                         m_positions.size());
}


// BoundarySystem ParticleSystem::createBoundary(const std::vector<Eigen::Vector3d> &positions,
//                                               double restDensity,
//                                               double particleMass,
//                                               double particleRadius)

// {
//     
// }





// double FluidSystem::getTimeCFL()
// {
//     auto sort_func = [](Eigen::Vector3d v1, Eigen::Vector3d v2)
//     {
//         return v1.norm() < v2.norm();
//     };

//     std::vector<Eigen::Vector3d>::iterator maxVelocity = std::max_element(velocities.begin(),
//                                                                           velocities.end(),
//                                                                           sort_func);
//     return 0.5 * (particleRadius / (*maxVelocity).norm());
// }


// void FluidSystem::setAccelerationsToGravity() {
//     auto grav_func = [](Eigen::Vector3d a)
//     {
//         return Eigen::Vector3d(0.0, -9.80665, 0.0);
//     };
//     std::transform(accelerations.begin(), accelerations.end(),
//                    accelerations.begin(), grav_func);
// }


// std::vector<Eigen::Vector3d> ParticleSystem::interpolatePositions(const std::vector<Eigen::Vector3d>& prevPositions,
//                                                                   const std::vector<Eigen::Vector3d>& curPositions,
//                                                                   const double prevTime, const double curTime,
//                                                                   const double desiredTime)
// {
//     double timeBtwnCurAndPrev = curTime - prevTime;
//     double timeBtwnDesiredAndPrev = desiredTime - prevTime;
//     double beta = timeBtwnDesiredAndPrev / timeBtwnCurAndPrev;
//     double alpha = 1.0 - beta;

//     auto inter_func =  [alpha, beta](Eigen::Vector3d prevPos, Eigen::Vector3d curPos)
//     {
//         return alpha * prevPos + beta * curPos;
//     };
    
//     std::vector<Eigen::Vector3d> interpolation(prevPositions);
//     std::transform(prevPositions.begin(), prevPositions.end(),
//                    curPositions.begin(), interpolation.begin(), inter_func);
    
//     return interpolation;
// }


// void FluidSystem::updateAccelerations(std::vector<BoundarySystem> boundaries,
//                                       const CompactNSearch::NeighborhoodSearch& nsearch)
// {
//     Eigen::Vector3d g = Eigen::Vector3d(0.0, -9.80665, 0.0);
//     for (size_t fpI = 0; fpI < this->positions.size(); fpI++) {
//         Eigen::Vector3d accelerationPressure = calculateAccelerationPressure(boundaries, nsearch, fpI);
//         Eigen::Vector3d accelerationViscosity = calculateAccelerationViscosity(boundaries, nsearch, fpI);

//         this->accelerations[fpI] = g + accelerationPressure + accelerationViscosity;
//     }
// }


// Eigen::Vector3d FluidSystem::calculateAccelerationPressure(std::vector<BoundarySystem> boundaries,
//                                                            const CompactNSearch::NeighborhoodSearch& nsearch,
//                                                            const size_t fpI)
// {
//     double particleDensityRatio = calculatePressure(this->densities[fpI], this->restDensity) / pow(this->densities[fpI], 2);
//     Eigen::Vector3d accelerationPressureFF = calculateAccelerationPressureFF(nsearch, fpI, particleDensityRatio);
//     Eigen::Vector3d accelerationPressureFS = calculateAccelerationPressureFS(boundaries, nsearch, fpI, particleDensityRatio);
//     return -accelerationPressureFF - accelerationPressureFS;
// }


// Eigen::Vector3d FluidSystem::calculateAccelerationPressureFF(const CompactNSearch::NeighborhoodSearch& nsearch,
//                                                              const size_t fpI, const double particleDensityRatio)
// {
//     const double h = Kernel::Parameter::TUNING * this->particleRadius * 2;
//     CompactNSearch::PointSet const& fluidPS = nsearch.point_set(this->id);
//     Eigen::Vector3d result = Eigen::Vector3d::Zero();
    
//     for (size_t fpN = 0; fpN < fluidPS.n_neighbors(this->id, fpI); fpN++) {
//         const unsigned int fnI = fluidPS.neighbor(this->id, fpI, fpN);

//         double neighborDensityRatio = calculatePressure(this->densities[fnI], this->restDensity) / pow(this->densities[fnI], 2);
//         double pressureTerm = particleDensityRatio + neighborDensityRatio;
//         result += this->particleMass * pressureTerm * Kernel::CubicSpline::gradWeight(this->positions[fpI], this->positions[fnI], h);
//     }
//     return result;
// }


// Eigen::Vector3d FluidSystem::calculateAccelerationPressureFS(std::vector<BoundarySystem> boundaries,
//                                                              const CompactNSearch::NeighborhoodSearch& nsearch,
//                                                              const size_t fpI, const double particleDensityRatio)
// {
//     Eigen::Vector3d result = Eigen::Vector3d::Zero();
//     if (boundaries.size() > 0) {
//         for (const BoundarySystem &boundary : boundaries) {
//             const double h = Kernel::Parameter::TUNING * this->particleRadius * 2;
//             CompactNSearch::PointSet const& fluidPS = nsearch.point_set(this->id);
//             for (size_t fpN = 0; fpN < fluidPS.n_neighbors(boundary.id, fpI); fpN++) {
//                 const unsigned int fnI = fluidPS.neighbor(boundary.id, fpI, fpN);
//                 result += this->restDensity * boundary.volumes[fnI]
//                     * particleDensityRatio
//                     * Kernel::CubicSpline::gradWeight(this->positions[fpI], boundary.positions[fnI], h);
//             }
//         }
//     }
//     return result;
// }


// Eigen::Vector3d FluidSystem::calculateAccelerationViscosity(std::vector<BoundarySystem> boundaries,
//                                                             const CompactNSearch::NeighborhoodSearch& nsearch,
//                                                             const size_t fpI)
// {
//     Eigen::Vector3d accelerationViscosityFF = calculateAccelerationViscosityFF(nsearch, fpI);
//     Eigen::Vector3d accelerationViscosityFS = calculateAccelerationViscosityFS(boundaries, nsearch, fpI);
    
//     return 2.0 * this->viscosity * accelerationViscosityFF + accelerationViscosityFS;
// }

// Eigen::Vector3d FluidSystem::calculateAccelerationViscosityFF(const CompactNSearch::NeighborhoodSearch& nsearch,
//                                                               const size_t fpI)
// {
//     const double h = Kernel::Parameter::TUNING * this->particleRadius * 2;
//     CompactNSearch::PointSet const& fluidPS = nsearch.point_set(this->id);
//     Eigen::Vector3d result = Eigen::Vector3d::Zero();
//     for (size_t fpN = 0; fpN < fluidPS.n_neighbors(this->id, fpI); fpN++) {
//         const unsigned int fnI = fluidPS.neighbor(this->id, fpI, fpN);

//         double neighborVol = this->particleMass / this->densities[fnI];
//         Eigen::Vector3d velocityDiff = this->velocities[fpI] - this->velocities[fnI];
//         Eigen::Vector3d posDiff = this->positions[fpI] - this->positions[fnI];
//         result += neighborVol * velocityDiff * posDiff.transpose()
//             * Kernel::CubicSpline::gradWeight(this->positions[fpI], this->positions[fnI], h)
//             * pow(pow(posDiff.norm(), 2) + 0.01 * pow(h, 2), -1);
//     }
//     return result;
// }


// Eigen::Vector3d FluidSystem::calculateAccelerationViscosityFS(std::vector<BoundarySystem> boundaries,
//                                                               const CompactNSearch::NeighborhoodSearch& nsearch,
//                                                               const size_t fpI)
// {
//     Eigen::Vector3d result = Eigen::Vector3d::Zero();
//     for (const BoundarySystem& boundary : boundaries) {
//         Eigen::Vector3d resultBoundary = Eigen::Vector3d::Zero();
//         const double h = Kernel::Parameter::TUNING * this->particleRadius * 2;
//         CompactNSearch::PointSet const& fluidPS = nsearch.point_set(this->id);
//         for (size_t fpN = 0; fpN < fluidPS.n_neighbors(boundary.id, fpI); fpN++) {
//             const unsigned int fnI = fluidPS.neighbor(boundary.id, fpI, fpN);

//             Eigen::Vector3d posDiff = this->positions[fpI] - boundary.positions[fnI];
//             resultBoundary += boundary.volumes[fnI]
//                 * this->velocities[fpI] * posDiff.transpose()
//                 * Kernel::CubicSpline::gradWeight(this->positions[fpI], boundary.positions[fnI], h)
//                 * pow(pow(posDiff.norm(), 2) + 0.01 * pow(h, 2), -1);
//         }
//         resultBoundary *= 2.0 * boundary.viscosity;
//         result += resultBoundary;
//     }
//     return result;
// }


// double ParticleSystem::calculatePressure(const double particleDensity,
//                                          const double restDensity)
// {
//     return std::max(particleDensity - restDensity, 0.0) * Parameter::B;
// }
