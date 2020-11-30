#pragma once
#include "particlesystem.h"
#include <CompactNSearch/CompactNSearch.h>

namespace learnSPH {
    class BoundarySystem;
    class FluidSystem : public ParticleSystem {
        friend class ParticleEmitter;
        
    private:
        FluidSystem() {}

    public:
        /* Estimate the densitiy of all fluid particles including the given
         * boundaries with the given neighbourhood information.
         * @param &boundaries - list of boundary systems.
         * @param &nsearch    - CompactNSearch neighbourhood information. */
        void estimateDensity(CompactNSearch::NeighborhoodSearch &nsearch,
                             const std::vector<BoundarySystem> &boundaries);
        void estimateDensity(CompactNSearch::NeighborhoodSearch &nsearch);


        // Setter & Getter
        std::vector<double> getDensities() const { return m_densities; }


    private:
        std::vector<double> m_densities; // Estimated particle densities
        std::vector<Eigen::Vector3d> m_velocities; // Particle Velocities
        std::vector<Eigen::Vector3d> m_accelerations; // Particle Accelerations

        
        
        /* /\* @return Calculates a timestep which satisfies the Courant Friedrich Levy (CFL) condition *\/ */
        /* double getTimeCFL(); */

        /* /\* Set the acceleration for every particle to g:= (0.0, -9.80665, 0.0) *\/ */
        /* void setAccelerationsToGravity(); */

        /* Eigen::Vector3d calculateAccelerationPressure(std::vector<BoundarySystem> boundaries, */
        /*                                               const CompactNSearch::NeighborhoodSearch& nsearch, */
        /*                                               const size_t fpI); */

        /* Eigen::Vector3d calculateAccelerationPressureFF(const CompactNSearch::NeighborhoodSearch& nsearch, */
        /*                                                 const size_t fpI, const double particleDensityRatio); */

        /* Eigen::Vector3d calculateAccelerationPressureFS(std::vector<BoundarySystem> boundaries, */
        /*                                                 const CompactNSearch::NeighborhoodSearch& nsearch, */
        /*                                                 const size_t fpI, const double particleDensityRatio); */

        /* Eigen::Vector3d calculateAccelerationViscosity(std::vector<BoundarySystem> boundaries, */
        /*                                                const CompactNSearch::NeighborhoodSearch& nsearch, */
        /*                                                const size_t fpI); */

        /* Eigen::Vector3d calculateAccelerationViscosityFF(const CompactNSearch::NeighborhoodSearch& nsearch, */
        /*                                                  const size_t fpI); */

        /* Eigen::Vector3d calculateAccelerationViscosityFS(std::vector<BoundarySystem> boundaries, */
        /*                                                  const CompactNSearch::NeighborhoodSearch& nsearch, */
        /*                                                  const size_t fpI); */
            
        /* void updateAccelerations(std::vector<BoundarySystem> boundaries, const CompactNSearch::NeighborhoodSearch& nsearch); */
    };
} // namespace learnSPH
