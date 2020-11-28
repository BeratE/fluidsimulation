#pragma once
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <CompactNSearch/CompactNSearch.h>

namespace learnSPH {
    namespace ParticleSystem {
        namespace Parameter {
            constexpr double B = 5000.0; // B form calculation of pressure for particles (see Assignment 2, Exercise 4)
        } //namespace Parameter
        struct ParticleSystem {
            long id = -1; // CompactNSearch PointSet ID
            double restDensity = 0.0;
            double particleMass = 0.0; // Constant particle density for all particles
            double particleRadius = 0.0; // (cubic root of volume) / 2
            double viscosity = 0.0; // Viscosity, should be different for fluid and boundary
            std::vector<Eigen::Vector3d> positions; // Particle Positions
        };

        struct BoundarySystem : ParticleSystem{
            std::vector<double> volumes; // Corrected particle volume
        };

        struct FluidSystem : ParticleSystem {
            std::vector<double> densities; // Estimated particle densities
            std::vector<Eigen::Vector3d> velocities; // Particle Velocities
            std::vector<Eigen::Vector3d> accelerations; // Particle Accelerations

            /*
             * @return Calculates a timestep which satisfies the Courant Friedrich Levy (CFL) condition 
             */
            double getTimeCFL();

            /*
             * Set the acceleration for every particle to g:= (0.0, -9.80665, 0.0)  
             */
            void setAccelerationsToGravity();

            Eigen::Vector3d calculateAccelerationPressure(std::vector<BoundarySystem> boundaries, const CompactNSearch::NeighborhoodSearch& nsearch, const size_t fpI);

            Eigen::Vector3d calculateAccelerationPressureFF(const CompactNSearch::NeighborhoodSearch& nsearch, const size_t fpI, const double particleDensityRatio);

            Eigen::Vector3d calculateAccelerationPressureFS(std::vector<BoundarySystem> boundaries, const CompactNSearch::NeighborhoodSearch& nsearch, const size_t fpI, const double particleDensityRatio);

            Eigen::Vector3d calculateAccelerationViscosity(std::vector<BoundarySystem> boundaries, const CompactNSearch::NeighborhoodSearch& nsearch, const size_t fpI);

            Eigen::Vector3d calculateAccelerationViscosityFF(const CompactNSearch::NeighborhoodSearch& nsearch, const size_t fpI);

            Eigen::Vector3d calculateAccelerationViscosityFS(std::vector<BoundarySystem> boundaries, const CompactNSearch::NeighborhoodSearch& nsearch, const size_t fpI);
            
            void updateAccelerations(std::vector<BoundarySystem> boundaries, const CompactNSearch::NeighborhoodSearch& nsearch);
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

        /* Sample positions in a triangle specified by 3d points a, b , c with sampling 
         * distance apart from each other (hexagonal). 
         * @param a, b , c         - 3d points of the triangle. 
         * @param samplingDistance - distance the sample points are apart. 
         * @return Vector of sampled 3d points within the triangle. */
        std::vector<Eigen::Vector3d> samplePositionsTriangle(Eigen::Vector3d a,
                                                             Eigen::Vector3d b,
                                                             Eigen::Vector3d c,
                                                             double samplingDistance);

        /* Sample a hollow box spanning from the bottom left to the top right corner.
         * @param bottomLeft - bottom left box corner.
         * @param topRight   - top right box corner diagonal to bottom left.
         * @param samplingDistance - particles are sampling distance apart.
         * @return vector of positions of sample points. */
        std::vector<Eigen::Vector3d> samplePositionsBox(Eigen::Vector3d bottomLeft,
                                                        Eigen::Vector3d topRight,
                                                        double samplingDistance);

        /* Construct a boundary from the given positions, rest density and mass */
        BoundarySystem createBoundary(const std::vector<Eigen::Vector3d> &positions,
                                      double restDensity, double particleMass,
                                      double particleRadius);

        /* Correct the volume of each particle in the given boundary system.  */
        void correctBoundaryVolume(BoundarySystem &boundary);

        /*
         * Estimates the positions of particles at desiredTime by interpolating the positions at prevTime and curTime
         * @param &prevPositions - Previous positions of the particles
         * @param &curPositions - Current positions of the particles
         * @param prevTime - Time at which the prevPositions were determined
         * @param curTime - Time at which the curPositions were determined
         * @param desiredTime - Time at which the positions of the particles should be estimated. Requires prevTime <= desiredTime <= curTime.
         * @returns interpolated positions of particles at desiredTime
         */
        std::vector<Eigen::Vector3d> interpolatePositions(const std::vector<Eigen::Vector3d>& prevPositions, const std::vector<Eigen::Vector3d>& curPositions, const double prevTime, const double curTime, const double desiredTime);
        
        double calculatePressure(const double particleDensity, const double restDensity);

    } // namespace ParticleSystem
} // namespace learnSPH
