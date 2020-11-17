#pragma once
#include <vector>
#include <Eigen/Dense>
#include <CompactNSearch/CompactNSearch.h>

namespace learnSPH {
    namespace ParticleSystem {

        struct ParticleSystem {
            unsigned int id;        // CompactNSearch PointSet ID
            double restDensity;
            std::vector<Eigen::Vector3d> positions; // Particle Positions
        };

        struct FluidSystem : ParticleSystem {
            double particleMass;      
            std::vector<double> densities; // Estimated particle densities
        };

        struct BoundarySystem : ParticleSystem{
            std::vector<double> mass;
        };

        /* Estimate the densitiy of all fluid particles including the given boundaries
         *  with the given neighbourhood information.
         * @param &fluid    - Fluid system struct. Densitiy is written back into the struct.
         * @param &boundary - Boundary system.
         * @param &nsearch  - CompactNSearch neighbourhood information */
        void estimateDensity(FluidSystem& fluid, BoundarySystem& boundary,
                             CompactNSearch::NeighborhoodSearch &nsearch);
        
    } // namespace ParticleSystem
} // namespace learnSPH
