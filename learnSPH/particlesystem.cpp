#include "particlesystem.h"
#include "kernel.h"

using namespace learnSPH;
using namespace learnSPH::ParticleSystem;

FluidSystem learnSPH::ParticleSystem::sampleBox(Eigen::Vector3d bottomLeft,
                                                Eigen::Vector3d topRight,
                                                double fluidRestDensity,
                                                double samplingDistance)
{
    const Eigen::Vector3d diagonal = topRight - bottomLeft;
    const Eigen::Vector3d diagonalSign = diagonal.array() / diagonal.array();
    const Eigen::Vector3d samplingDir = samplingDistance * diagonalSign;
    const Eigen::Vector3d numSamples = Eigen::floor(Eigen::abs(diagonal.array())/samplingDistance);
    
    FluidSystem fluidParticles;
    fluidParticles.restDensity = fluidRestDensity;
    fluidParticles.particleMass = pow(samplingDistance, 3) * fluidRestDensity;

    for (size_t x = 0; x < (size_t)numSamples[0]; x++) {
        for (size_t y = 0; y < (size_t)numSamples[1]; y++) {
            for (size_t z = 0; z < (size_t)numSamples[2]; z++) {
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

void learnSPH::ParticleSystem::estimateDensity(FluidSystem &fluid,
                                               CompactNSearch::NeighborhoodSearch &nsearch)
{
    std::vector<BoundarySystem> emptyList;
    estimateDensity(fluid, emptyList, nsearch);
}

void learnSPH::ParticleSystem::estimateDensity(FluidSystem& fluid,
                                               std::vector<BoundarySystem>& boundaries,
                                               CompactNSearch::NeighborhoodSearch& nsearch)
{
    const double h = Kernel::Parameter::TUNING * fluid.particleRadius * 2;
    
    CompactNSearch::PointSet const& fluidPS = nsearch.point_set(fluid.id);
    fluid.densities.resize(fluidPS.n_points());
    
    for (size_t fpI = 0; fpI < fluidPS.n_points(); fpI++) {
        const Eigen::Vector3d &fpPos = fluid.positions[fpI];
        double particleDensity = 0.0;

        // Fluid contribution
        /* Note: The current fluid particle itself is part of its own neighborhood. */
        particleDensity += fluid.particleMass * fluidPS.n_neighbors(fluid.id, fpI +1);
        particleDensity += Kernel::CubicSpline::weight(fpPos, fpPos, h);
        for (size_t fpN = 0; fpN < fluidPS.n_neighbors(fluid.id, fpI); fpN++) {
            const unsigned int fnI = fluidPS.neighbor(fluid.id, fpI, fpN);
            particleDensity += Kernel::CubicSpline::weight(fpPos, fluid.positions[fnI], h);
        }

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
