#include "particlesystem.h"
#include "kernel.h"

using namespace learnSPH;
using namespace learnSPH::ParticleSystem;

void estimateDensity(FluidSystem& fluid, std::vector<BoundarySystem>& boundaries,
                     CompactNSearch::NeighborhoodSearch& nsearch)
{
    const double h = 0.2;

   
    CompactNSearch::PointSet const& pointSet = nsearch.point_set(fluid.id);
    fluid.densities.resize(pointSet.n_points());
    
    for (size_t i = 0; i < pointSet.n_points(); i++) {
        double density = 0.0;
        Eigen::Vector3d &particlePos = fluid.positions[i];
        
        // Fluid contribution
        for (size_t j = 0; j < pointSet.n_neighbors(fluid.id, i); j++) {
            const unsigned int nid = pointSet.neighbor(fluid.id, i, j);
            density += fluid.particleMass
                * Kernel::CubicSpline::weight(particlePos, fluid.positions[nid], h);
        }

        // Boundary contribution
        for (const BoundarySystem &boundary : boundaries) {
            for (size_t j = 0; j < pointSet.n_neighbors(boundary.id, i); j++) {
                const unsigned int nid = pointSet.neighbor(boundary.id, i, j);
                density += boundary.mass[nid]
                    * Kernel::CubicSpline::weight(particlePos, boundary.positions[nid], h);
                
            }
        }

        fluid.densities[i] = density;
    }
    
}
