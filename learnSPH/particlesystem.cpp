#include <functional>
#include "particlesystem.h"
#include "kernel.h"

using namespace learnSPH;
using namespace learnSPH::ParticleSystem;


FluidSystem
learnSPH::ParticleSystem::sampleFluidBox(Eigen::Vector3d bottomLeft,
                                         Eigen::Vector3d topRight,
                                         double fluidRestDensity,
                                         double samplingDistance)
{
    const Eigen::Vector3d diagonal = topRight - bottomLeft;
    const Eigen::Vector3d diagonalSign = diagonal.array() / diagonal.array().abs();
    const Eigen::Vector3d samplingDir = samplingDistance * diagonalSign;
    const Eigen::Vector3d numSamples = Eigen::floor(Eigen::abs(diagonal.array())/samplingDistance);
    
    FluidSystem fluidParticles;
    fluidParticles.restDensity = fluidRestDensity;
    fluidParticles.particleRadius = samplingDistance/2;
    fluidParticles.particleMass = pow(samplingDistance, 3) * fluidRestDensity;

    for (size_t x = 0; x < (size_t)numSamples[0]; x++) {
        for (size_t y = 0; y < (size_t)numSamples[1]; y++) {
            for (size_t z = 0; z < (size_t)numSamples[2]; z++) {
                // NOTE(DENNIS) - Kommen die letzten Samples hier nicht manchmal auch von au�erhalb des Volumens?
                Eigen::Vector3d samplingPos ((x + 0.5) * samplingDir[0],
                                             (y + 0.5) * samplingDir[1],
                                             (z + 0.5) * samplingDir[2]);
                samplingPos += bottomLeft;

                fluidParticles.positions.push_back(samplingPos);
            }
        }
    }

    return fluidParticles;
}


std::vector<Eigen::Vector3d>
learnSPH::ParticleSystem::samplePositionsTriangle(Eigen::Vector3d a,
                                                  Eigen::Vector3d b,
                                                  Eigen::Vector3d c,
                                                  double samplingDistance)
{   
    // Construct orthonormal base of the 2d plane with origin at a spanned by the triangle
    Eigen::Vector3d u = (b-a).normalized();
    Eigen::Vector3d n = u.cross((c-a).normalized()).normalized();
    Eigen::Vector3d v = n.cross(u).normalized();
    Eigen::Matrix<double, 3, 4> planeBasis;
    planeBasis.row(0) = Eigen::Vector4d(u(0), u(1), u(2), -u.dot(a));
    planeBasis.row(1) = Eigen::Vector4d(v(0), v(1), v(2), -v.dot(a));
    planeBasis.row(2) = Eigen::Vector4d(0.0, 0.0, 0.0, 1.0);
    
    // Convert triangle points to new 2d basis (discard homogenous coordinate)
    Eigen::Vector2d x = (planeBasis * a.homogeneous()).segment(0, 2);
    Eigen::Vector2d y = (planeBasis * b.homogeneous()).segment(0, 2);
    Eigen::Vector2d z = (planeBasis * c.homogeneous()).segment(0, 2);

    // Return the hessian line function of the edge from point a to point b,
    // offset by half the sampling distance from the opposite triangle point o.
    auto constructEdge = [samplingDistance](Eigen::Vector2d a,
                                            Eigen::Vector2d b,
                                            Eigen::Vector2d o) -> std::function<double(Eigen::Vector2d)> {
        Eigen::Vector2d normal = (b-a);
        double tmp = normal(0);
        normal(0) = normal(1);
        normal(1) = -tmp;
        normal.normalize();
        Eigen::Vector2d offset = ((a + (b - a)/2) - o).normalized() * samplingDistance/2;
        double distance = normal.dot(a + offset);
        return [=](Eigen::Vector2d point) mutable {
            return normal.dot(point) - distance;
        };
    };

    auto edgeXY = constructEdge(x, y, z);
    auto edgeYZ = constructEdge(y, z, x);
    auto edgeZX = constructEdge(z, x, y);
    
    // Generate points in the 2d rectancle enclosing the triangle with the
    // pineda algorithm
    Eigen::Vector2d rectMin = -(samplingDistance/2.0)*Eigen::Vector2d(1.0, 1.0); // rectangle (minX, minY)
    Eigen::Vector2d rectMax =  (samplingDistance/2.0)*Eigen::Vector2d(1.0, 1.0); // rectangle (maxX, maxY)
    rectMin(0) += (x(0) < y(0)) ? (x(0) < z(0) ? x(0) : z(0)) : (y(0) < z(0) ? y(0) : z(0));
    rectMin(1) += (x(1) < y(1)) ? (x(1) < z(1) ? x(1) : z(1)) : (y(1) < z(1) ? y(1) : z(1));
    rectMax(0) += (x(0) > y(0)) ? (x(0) > z(0) ? x(0) : z(0)) : (y(0) > z(0) ? y(0) : z(0));
    rectMax(1) += (x(1) > y(1)) ? (x(1) > z(1) ? x(1) : z(1)) : (y(1) > z(1) ? y(1) : z(1));

    // Sample points in a raster of integer steps of sampling distance
    std::vector<Eigen::Vector2d> samples;
    for (size_t i = 0; i < std::floor((rectMax(0) - rectMin(0)) / samplingDistance); i++) {
        for (size_t j = 0; j < std::floor((rectMax(1) - rectMin(1)) / samplingDistance); j++) {
            Eigen::Vector2d sample(rectMin(0) + i * samplingDistance,
                                   rectMin(1) + j * samplingDistance);
            if (edgeXY(sample) <= 0 && edgeYZ(sample) <= 0 && edgeZX(sample) <= 0)
                samples.push_back(sample);
        }
    }

    std::vector<Eigen::Vector3d> positions;
    // Transform sample points back to 3d standard basis
    for (const auto &s : samples) {
        Eigen::Vector3d sample_point = (planeBasis.transpose() * s.homogeneous()).segment(0, 3) + a;
        positions.push_back(sample_point);
    }

    return positions;
}

std::vector<Eigen::Vector3d>
learnSPH::ParticleSystem::samplePositionsBox(Eigen::Vector3d bottomLeft,
                                             Eigen::Vector3d topRight,
                                             double samplingDistance)
{
    std::vector<Eigen::Vector3d> positions;
    // Generate corner vertices
    Eigen::Vector3d dir = topRight - bottomLeft;
    std::vector<Eigen::Vector3d> vertices;
    vertices.push_back(Eigen::Vector3d(bottomLeft(0),        bottomLeft(1),        bottomLeft(2)));
    vertices.push_back(Eigen::Vector3d(bottomLeft(0),        bottomLeft(1)+dir(1), bottomLeft(2)));
    vertices.push_back(Eigen::Vector3d(bottomLeft(0)+dir(0), bottomLeft(1),        bottomLeft(2)));
    vertices.push_back(Eigen::Vector3d(bottomLeft(0)+dir(0), bottomLeft(1)+dir(1), bottomLeft(2)));
    vertices.push_back(Eigen::Vector3d(bottomLeft(0),        bottomLeft(1),        bottomLeft(2)+dir(1)));
    vertices.push_back(Eigen::Vector3d(bottomLeft(0),        bottomLeft(1)+dir(1), bottomLeft(2)+dir(1)));
    vertices.push_back(Eigen::Vector3d(bottomLeft(0)+dir(0), bottomLeft(1),        bottomLeft(2)+dir(1)));
    vertices.push_back(Eigen::Vector3d(bottomLeft(0)+dir(0), bottomLeft(1)+dir(1), bottomLeft(2)+dir(1)));

    std::vector<std::vector<Eigen::Vector3d>> v;
    v.push_back(samplePositionsTriangle(vertices[2], vertices[1], vertices[3], samplingDistance));
    v.push_back(samplePositionsTriangle(vertices[0], vertices[1], vertices[2], samplingDistance));
    v.push_back(samplePositionsTriangle(vertices[0], vertices[4], vertices[2], samplingDistance));
    v.push_back(samplePositionsTriangle(vertices[2], vertices[4], vertices[6], samplingDistance));
    v.push_back(samplePositionsTriangle(vertices[0], vertices[4], vertices[1], samplingDistance));
    v.push_back(samplePositionsTriangle(vertices[1], vertices[4], vertices[5], samplingDistance));
    v.push_back(samplePositionsTriangle(vertices[1], vertices[5], vertices[3], samplingDistance));
    v.push_back(samplePositionsTriangle(vertices[3], vertices[5], vertices[7], samplingDistance));
    v.push_back(samplePositionsTriangle(vertices[2], vertices[6], vertices[3], samplingDistance));
    v.push_back(samplePositionsTriangle(vertices[3], vertices[6], vertices[7], samplingDistance));
    v.push_back(samplePositionsTriangle(vertices[4], vertices[5], vertices[6], samplingDistance));
    v.push_back(samplePositionsTriangle(vertices[6], vertices[5], vertices[7], samplingDistance));
    
    for (const auto &p : v) {
        positions.insert(positions.end(), p.begin(), p.end());
    }

    return positions;
}

void
learnSPH::ParticleSystem::estimateFluidDensity(FluidSystem &fluid,
                                               const CompactNSearch::NeighborhoodSearch &nsearch)
{
    estimateFluidDensity(fluid, std::vector<BoundarySystem>(), nsearch);
}


void
learnSPH::ParticleSystem::estimateFluidDensity(FluidSystem& fluid,
                                               const std::vector<BoundarySystem>& boundaries,
                                               const CompactNSearch::NeighborhoodSearch& nsearch)
{
    // smoothing length
    const double h = Kernel::Parameter::TUNING * fluid.particleRadius * 2;
    // get neighborhood information of fluid particle point set
    CompactNSearch::PointSet const& fluidPS = nsearch.point_set(fluid.id);
    fluid.densities.resize(fluidPS.n_points());
    // iterate fluid particles
    for (size_t fpI = 0; fpI < fluidPS.n_points(); fpI++) {
        const Eigen::Vector3d &fpPos = fluid.positions[fpI];
        double particleDensity = 0.0;

        // Fluid contribution
        /* Note: The current fluid particle itself is part of its own neighborhood. */
        particleDensity += Kernel::CubicSpline::weight(fpPos, fpPos, h);
        for (size_t fpN = 0; fpN < fluidPS.n_neighbors(fluid.id, fpI); fpN++) {
            const unsigned int fnI = fluidPS.neighbor(fluid.id, fpI, fpN);
            particleDensity += Kernel::CubicSpline::weight(fpPos, fluid.positions[fnI], h);
        }
        particleDensity *= fluid.particleMass;

        // Boundary contribution
        for (const BoundarySystem &boundary : boundaries) {
            for (size_t bpI = 0; bpI < fluidPS.n_neighbors(boundary.id, fpI); bpI++) {
                const unsigned int bnI = fluidPS.neighbor(boundary.id, fpI, bpI);
                particleDensity += boundary.particleMasses[bnI]
                    * Kernel::CubicSpline::weight(fpPos, boundary.positions[bnI], h);
            }
        }

        fluid.densities[fpI] = particleDensity;
    }  
}

double FluidSystem::getTimeCFL() {
    std::vector<Eigen::Vector3d>::iterator maxVelocity = std::max_element(velocities.begin(), velocities.end(), [](Eigen::Vector3d v1, Eigen::Vector3d v2) {
        return v1.norm() < v2.norm(); });
    return 0.5 * (particleRadius / (*maxVelocity).norm());
}
