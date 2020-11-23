#pragma once
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <CompactNSearch/CompactNSearch.h>

namespace learnSPH {
    namespace ParticleSystem {

        struct ParticleSystem {
            long id = -1; // CompactNSearch PointSet ID
            double restDensity = 0.0;
            double particleRadius = 0.0; // (cubic root of volume) / 2
            std::vector<Eigen::Vector3d> positions; // Particle Positions
        };

        struct FluidSystem : ParticleSystem {
            double particleMass = 0.0; // Constant particle density for all fluid particles
            std::vector<double> densities; // Estimated particle densities
        };

        struct BoundarySystem : ParticleSystem{
            std::vector<double> particleMasses; // Each particle has its own mass
        };


        /* Uniformly sample an axis box with fluid particles. 
         * @param bottomLeft       - bottom left corner of 3d axis aligned box eg. (0,0,0). 
         * @param topRight         - top right corner diagonal to bottom left eg. (1,1,1). 
         * @param fluidRestDensity - rest density of particle system. 
         * @param samplingDistance - distance between particle samples. 
         *                           Effective particle diameter.
         * @return FuidSystem particles in the axis aligned box. */        
        FluidSystem sampleFluidBox(Eigen::Vector3d bottomLeft, Eigen::Vector3d topRight,
                                   double fluidRestDensity, double samplingDistance);

        /* Sample positions in a triangle specified by 3d points a, b , c with sampling 
         * distance apart from each other. 
         * @param a, b , c         - 3d points of the triangle. 
         * @param samplingDistance - distance the sample points are apart. 
         * @return Vector of sampled 3d points within the triangle. */
        std::vector<Eigen::Vector3d> samplePositionsTriangle(Eigen::Vector3d a,
                                                             Eigen::Vector3d b,
                                                             Eigen::Vector3d c,
                                                             double samplingDistance);
        
        std::vector<Eigen::Vector3d> samplePositionsBox(Eigen::Vector3d bottomLeft,
                                                        Eigen::Vector3d topRight,
                                                        double samplingDistance);
        
        /* Estimate the densitiy of all fluid particles including the given boundaries
         *  with the given neighbourhood information.
         * @param &fluid      - Fluid system struct. Densitiy is written back into the struct.
         * @param &boundaries - list of boundary systems.
         * @param &nsearch    - CompactNSearch neighbourhood information. */
        void estimateFluidDensity(FluidSystem& fluid,
                             const std::vector<BoundarySystem>& boundaries,
                             const CompactNSearch::NeighborhoodSearch &nsearch);
        void estimateFluidDensity(FluidSystem& fluid,
                             const CompactNSearch::NeighborhoodSearch &nsearch);
        
    } // namespace ParticleSystem
} // namespace learnSPH
