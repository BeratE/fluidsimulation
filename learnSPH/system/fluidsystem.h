#pragma once
#include "particlesystem.h"
#include <CompactNSearch/CompactNSearch.h>

namespace learnSPH::System {
    class BoundarySystem;
    class FluidSystem : public ParticleSystem {
        friend class ParticleEmitter;

    public:
        FluidSystem(double radius, double density, size_t size, bool fill = true);

        void updatePressures(double stiffness);
        void updateDensities(const std::vector<BoundarySystem> &boundaries);
        void updateAccelerations(const std::vector<BoundarySystem> &boundaries,
                                 bool pressure = true, bool viscosity = true, bool external = true); 
        
        // Setter & Getter
        double getParticleDensity(size_t i) const { return m_densities[i]; }
        double getParticlePressure(size_t i) const { return m_pressures[i]; }
        Eigen::Vector3d getParticleAcc(size_t i) const {return m_accelerations[i];}

        const std::vector<double> &getDensities() const { return m_densities; }
        const std::vector<double> &getPressures() const { return m_pressures; }
        const std::vector<Vector3d> &getAccelerations() const { return m_accelerations; }
        const std::vector<double> &getNormalizedDensities() const { return m_normalizedDensities; }

      private:        
        Vector3d particlePressureAcc(size_t i, const std::vector<BoundarySystem> &boundaries);
        Vector3d particleViscosityAcc(size_t i, const std::vector<BoundarySystem> &boundaries);       
        
        std::vector<double> m_pressures; // last updated particle pressures
        std::vector<double> m_densities; // last updated particle densities
        std::vector<Vector3d> m_accelerations; // last updated accelerations
    };
} // namespace learnSPH::System
