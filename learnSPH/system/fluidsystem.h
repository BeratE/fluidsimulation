#pragma once
#include "particlesystem.h"
#include <CompactNSearch/CompactNSearch.h>

namespace learnSPH::System {
    class BoundarySystem;
    class FluidSystem : public ParticleSystem {
        friend class ParticleEmitter;

    public:        
        FluidSystem() {}
        FluidSystem(size_t size, bool fill = true);
        FluidSystem(double particleRadius, size_t size, bool fill = true);

        void initKernelLookupTable();
        
        void updateDensities(CompactNSearch::NeighborhoodSearch &nsearch,
                             const std::vector<BoundarySystem> &boundaries);
        void updateDensities(CompactNSearch::NeighborhoodSearch &nsearch);
        void updatePressures(double stiffness);
        void updateAccelerations(CompactNSearch::NeighborhoodSearch& nsearch,
                                 const std::vector<BoundarySystem> &boundaries);
        
        // Setter & Getter
        double getParticleDensity(size_t i) {return m_densities[i]; }
        double getParticlePressure(size_t i) {return m_pressures[i]; }

        const std::vector<double>& getDensities() const { return m_densities; }
        const std::vector<double>& getPressires() const { return m_pressures; }


    private:        
        Eigen::Vector3d particlePressureAcc(size_t index,
                                            CompactNSearch::NeighborhoodSearch& nsearch,
                                            const std::vector<BoundarySystem> &boundaries);

        Eigen::Vector3d particleViscosityAcc(size_t index,
                                             CompactNSearch::NeighborhoodSearch &nsearch,
                                             const std::vector<BoundarySystem> &boundaries);
                   

        
        std::vector<double> m_pressures;
        std::vector<double> m_densities;

        bool m_isLookupTableInit = false;
        Kernel::CubicSpline::Table m_kernelLookup;
    };
} // namespace learnSPH::System
