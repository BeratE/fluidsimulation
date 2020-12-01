#pragma once
#include "particlesystem.h"
#include "kernel.h"
#include <CompactNSearch/CompactNSearch.h>

namespace learnSPH {
    class BoundarySystem;
    class FluidSystem : public ParticleSystem {
        friend class ParticleEmitter;

    public:
        
        FluidSystem() {}
        FluidSystem(size_t size, bool fill = true);

        void initTable();
        
        void estimateDensity(CompactNSearch::NeighborhoodSearch &nsearch,
                             const std::vector<BoundarySystem> &boundaries);
        void estimateDensity(CompactNSearch::NeighborhoodSearch &nsearch);
        void updatePressure();
        void updateAcceleration(CompactNSearch::NeighborhoodSearch& nsearch,
                                 const std::vector<BoundarySystem> &boundaries);
        
        // Setter & Getter
        void setStiffness(double value) { m_stiffness = value; }
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
                   

        double m_stiffness = 1000;
        std::vector<double> m_pressures;
        std::vector<double> m_densities; // Estimated particle densities

        Kernel::CubicSpline::Table m_kernelTable;
    };
} // namespace learnSPH
