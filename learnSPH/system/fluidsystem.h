#pragma once
#include "particlesystem.h"
#include <CompactNSearch/CompactNSearch.h>

namespace learnSPH::System {
    class BoundarySystem;
    class FluidSystem : public ParticleSystem {
        friend class ParticleEmitter;

    public:
        FluidSystem(double radius, double density, size_t size, bool fill = true);

        void clearTensionForces();
        void clearAdhesionForces();

        void addParticleTensionForce(size_t i, Vector3d force);
        void addParticleAdhesionForce(size_t i, Vector3d force);

        void updatePressures(double stiffness);
        void updateDensities(const std::vector<BoundarySystem> &boundaries);
        void updateAccelerations(const std::vector<BoundarySystem> &boundaries,
                                 bool pressure = true, bool viscosity = true, bool external = true, bool tension = true, bool adhesion = true);
        void updateNormals(const double c);
        
        // Setter & Getter
        double getParticleDensity(size_t i) const { return m_densities[i]; }
        double getParticlePressure(size_t i) const { return m_pressures[i]; }
        Eigen::Vector3d getParticleAcc(size_t i) const {return m_accelerations[i];}
        Eigen::Vector3d getParticleNormal(size_t i) const { return m_normals[i]; }

        const std::vector<double> &getDensities() const { return m_densities; }
        const std::vector<double> &getPressures() const { return m_pressures; }
        const std::vector<Vector3d> &getAccelerations() const { return m_accelerations; }
        const std::vector<double> &getNormalizedDensities() const { return m_normalizedDensities; }
        const std::vector<Eigen::Vector3d>& getNormals() const { return m_normals; }
        const std::vector<Vector3d>& getTensionForces() const { return m_tensionForces; }
        const std::vector<Vector3d>& getAdhesionForces() const { return m_adhesionForces; }

      private:        
        Vector3d particlePressureAcc(size_t i, const std::vector<BoundarySystem> &boundaries);
        Vector3d particleViscosityAcc(size_t i, const std::vector<BoundarySystem> &boundaries);       
        Eigen::Vector3d normal(const size_t i, const double c);
        
        std::vector<double> m_pressures; // last updated particle pressures
        std::vector<double> m_densities; // last updated particle densities
        std::vector<Vector3d> m_accelerations; // last updated accelerations
        std::vector<Vector3d> m_normals; // normals for surface tension calculations
        std::vector<Vector3d> m_tensionForces; // forces that result from surface tension
        std::vector<Vector3d> m_adhesionForces; // forces that result form adhesion to boundaries
    };
} // namespace learnSPH::System
