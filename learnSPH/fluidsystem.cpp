#include "fluidsystem.h"
#include "boundarysystem.h"
#include "kernel.h"

using namespace learnSPH;
using namespace learnSPH::Kernel;
using namespace CompactNSearch;


void FluidSystem::estimateDensity(NeighborhoodSearch &nsearch)
{
    estimateDensity(nsearch, std::vector<BoundarySystem>());
}

void FluidSystem::estimateDensity(NeighborhoodSearch &nsearch,
                                  const std::vector<BoundarySystem> &boundaries)
{
    // get neighborhood information of fluid particle point set
    CompactNSearch::PointSet const& fluidPS = nsearch.point_set(m_pointSetID);
    m_densities.resize(fluidPS.n_points());
    
    // iterate fluid particles
    for (size_t fpI = 0; fpI < fluidPS.n_points(); fpI++) {
        const Eigen::Vector3d &fpPos = m_positions[fpI];

        // Fluid contribution
        /* The current fluid particle itself is part of its own neighborhood. */
        double fluidDensity = 0.0;
        fluidDensity += CubicSpline::weight(fpPos, fpPos,
                                            smoothingLength());
        for (size_t fpN = 0; fpN < fluidPS.n_neighbors(m_pointSetID, fpI); fpN++) {
            const unsigned int fnI = fluidPS.neighbor(m_pointSetID, fpI, fpN);
            fluidDensity += CubicSpline::weight(fpPos, m_positions[fnI],
                                                smoothingLength());
        }
        fluidDensity *= particleMass();

        // Boundary contribution
        double boundaryDensity = 0.0;
        for (const BoundarySystem &boundary : boundaries) {
            double density = 0.0;
            size_t n_neighbors = fluidPS.n_neighbors(boundary.getPointSetID(), fpI);
            for (size_t bpI = 0; bpI < n_neighbors; bpI++) {
                const unsigned int bnI = fluidPS.neighbor(boundary.getPointSetID(),
                                                          fpI, bpI);
                density += boundary.getParticleVolume(bnI)
                    * CubicSpline::weight(fpPos, boundary.getParticlePos(bnI),
                                          smoothingLength());
            }
            density *= boundary.getRestDensity();
            boundaryDensity += density;
        }

        m_densities[fpI] = fluidDensity + boundaryDensity;
    }  
}

