#include <functional>
#include "particlesystem.h"
#include "kernel.h"

using namespace learnSPH;

double ParticleSystem::smoothingLength()
{
    return 2 * m_particleRadius * Kernel::Parameter::TUNING;
}






// std::vector<Eigen::Vector3d> ParticleSystem::samplePositionsTriangle(Eigen::Vector3d a,
//                                                                      Eigen::Vector3d b,
//                                                                      Eigen::Vector3d c,
//                                                                      double samplingDistance)
// {   
//     // Construct orthonormal base of the 2d plane with origin at a spanned by the triangle
//     Eigen::Vector3d u = (b-a).normalized();
//     Eigen::Vector3d n = u.cross((c-a).normalized()).normalized();
//     Eigen::Vector3d v = n.cross(u).normalized();
//     Eigen::Matrix<double, 3, 4> planeBasis;
//     planeBasis.row(0) = Eigen::Vector4d(u(0), u(1), u(2), -u.dot(a));
//     planeBasis.row(1) = Eigen::Vector4d(v(0), v(1), v(2), -v.dot(a));
//     planeBasis.row(2) = Eigen::Vector4d(0.0, 0.0, 0.0, 1.0);
    
//     // Convert triangle points to new 2d basis (discard homogenous coordinate)
//     Eigen::Vector2d x = (planeBasis * a.homogeneous()).segment(0, 2);
//     Eigen::Vector2d y = (planeBasis * b.homogeneous()).segment(0, 2);
//     Eigen::Vector2d z = (planeBasis * c.homogeneous()).segment(0, 2);

//     // Return the hessian line function of the edge from point a to point b,
//     // offset by half the sampling distance from the opposite triangle point o.
//     auto constructEdge = [samplingDistance](Eigen::Vector2d a,
//                                             Eigen::Vector2d b,
//                                             Eigen::Vector2d o)
//         -> std::function<double(Eigen::Vector2d)> {
//         Eigen::Vector2d normal = (b-a);
//         double tmp = normal(0);
//         normal(0) = normal(1);
//         normal(1) = -tmp;
//         normal.normalize();
//         Eigen::Vector2d offset = ((a + (b - a)/2) - o).normalized() * samplingDistance/2;
//         double distance = normal.dot(a + offset);
//         return [=](Eigen::Vector2d point) mutable {
//             return normal.dot(point) - distance;
//         };
//     };

//     auto edgeXY = constructEdge(x, y, z);
//     auto edgeYZ = constructEdge(y, z, x);
//     auto edgeZX = constructEdge(z, x, y);
    
//     // Generate points in the 2d rectancle enclosing the triangle with the
//     // pineda algorithm
//     Eigen::Vector2d rectMin = -(samplingDistance/2.0)*Eigen::Vector2d(1.0, 1.0); // rectangle (minX, minY)
//     Eigen::Vector2d rectMax =  (samplingDistance/2.0)*Eigen::Vector2d(1.0, 1.0); // rectangle (maxX, maxY)
//     rectMin(0) += (x(0) < y(0)) ? (x(0) < z(0) ? x(0) : z(0)) : (y(0) < z(0) ? y(0) : z(0));
//     rectMin(1) += (x(1) < y(1)) ? (x(1) < z(1) ? x(1) : z(1)) : (y(1) < z(1) ? y(1) : z(1));
//     rectMax(0) += (x(0) > y(0)) ? (x(0) > z(0) ? x(0) : z(0)) : (y(0) > z(0) ? y(0) : z(0));
//     rectMax(1) += (x(1) > y(1)) ? (x(1) > z(1) ? x(1) : z(1)) : (y(1) > z(1) ? y(1) : z(1));

//     // Sample points in a hexagonal raster of
//     // integer steps of sampling distance in j direction and sampling distance/sqrt(2) in i direction
//     std::vector<Eigen::Vector2d> samples;
//     for (size_t i = 0; i < std::floor((rectMax(0) - rectMin(0)) / (samplingDistance/sqrt(2))); i++) {
//         for (size_t j = 0; j < std::floor((rectMax(1) - rectMin(1)) / samplingDistance); j++) {
//             Eigen::Vector2d sample(rectMin(0) + i * (samplingDistance/sqrt(2)),
//                                    rectMin(1) + j * samplingDistance);
//             sample(1) += ((i%2) * samplingDistance*sqrt(2));
//             if (edgeXY(sample) <= 0 && edgeYZ(sample) <= 0 && edgeZX(sample) <= 0)
//                 samples.push_back(sample);
//         }
//     }

//     std::vector<Eigen::Vector3d> positions;
//     // Transform sample points back to 3d standard basis
//     for (const auto &s : samples) {
//         Eigen::Vector3d sample_point = (planeBasis.transpose() * s.homogeneous()).segment(0, 3) + a;
//         positions.push_back(sample_point);
//     }

//     return positions;
// }


// std::vector<Eigen::Vector3d> ParticleSystem::samplePositionsBox(Eigen::Vector3d bottomLeft,
//                                                                 Eigen::Vector3d topRight,
//                                                                 double samplingDistance)
// {
//     std::vector<Eigen::Vector3d> positions;
//     // Generate corner vertices
//     Eigen::Vector3d dir = topRight - bottomLeft;
//     std::vector<Eigen::Vector3d> vertices;
//     vertices.push_back(Eigen::Vector3d(bottomLeft(0),        bottomLeft(1),        bottomLeft(2)));
//     vertices.push_back(Eigen::Vector3d(bottomLeft(0),        bottomLeft(1)+dir(1), bottomLeft(2)));
//     vertices.push_back(Eigen::Vector3d(bottomLeft(0)+dir(0), bottomLeft(1),        bottomLeft(2)));
//     vertices.push_back(Eigen::Vector3d(bottomLeft(0)+dir(0), bottomLeft(1)+dir(1), bottomLeft(2)));
//     vertices.push_back(Eigen::Vector3d(bottomLeft(0),        bottomLeft(1),        bottomLeft(2)+dir(1)));
//     vertices.push_back(Eigen::Vector3d(bottomLeft(0),        bottomLeft(1)+dir(1), bottomLeft(2)+dir(1)));
//     vertices.push_back(Eigen::Vector3d(bottomLeft(0)+dir(0), bottomLeft(1),        bottomLeft(2)+dir(1)));
//     vertices.push_back(Eigen::Vector3d(bottomLeft(0)+dir(0), bottomLeft(1)+dir(1), bottomLeft(2)+dir(1)));

//     std::vector<std::vector<Eigen::Vector3d>> v;
//     v.push_back(samplePositionsTriangle(vertices[2], vertices[1], vertices[3], samplingDistance));
//     v.push_back(samplePositionsTriangle(vertices[0], vertices[1], vertices[2], samplingDistance));
//     v.push_back(samplePositionsTriangle(vertices[0], vertices[4], vertices[2], samplingDistance));
//     v.push_back(samplePositionsTriangle(vertices[2], vertices[4], vertices[6], samplingDistance));
//     v.push_back(samplePositionsTriangle(vertices[0], vertices[4], vertices[1], samplingDistance));
//     v.push_back(samplePositionsTriangle(vertices[1], vertices[4], vertices[5], samplingDistance));
//     v.push_back(samplePositionsTriangle(vertices[1], vertices[5], vertices[3], samplingDistance));
//     v.push_back(samplePositionsTriangle(vertices[3], vertices[5], vertices[7], samplingDistance));
//     v.push_back(samplePositionsTriangle(vertices[2], vertices[6], vertices[3], samplingDistance));
//     v.push_back(samplePositionsTriangle(vertices[3], vertices[6], vertices[7], samplingDistance));
//     v.push_back(samplePositionsTriangle(vertices[4], vertices[5], vertices[6], samplingDistance));
//     v.push_back(samplePositionsTriangle(vertices[6], vertices[5], vertices[7], samplingDistance));
    
//     for (const auto &p : v) {
//         positions.insert(positions.end(), p.begin(), p.end());
//     }

//     return positions;
// }


// BoundarySystem ParticleSystem::createBoundary(const std::vector<Eigen::Vector3d> &positions,
//                                               double restDensity,
//                                               double particleMass,
//                                               double particleRadius)

// {
//     BoundarySystem boundary;
//     boundary.restDensity = restDensity;
//     boundary.particleMass = particleMass;
//     boundary.particleRadius = particleRadius;
//     boundary.positions.insert(boundary.positions.end(), positions.begin(), positions.end());
//     correctBoundaryVolume(boundary);
//     return boundary;
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
