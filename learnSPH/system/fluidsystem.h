#pragma once
#include "particlesystem.h"
#include <CompactNSearch/CompactNSearch.h>

namespace learnSPH::System {
    class BoundarySystem;
    class FluidSystem : public ParticleSystem {
        friend class ParticleEmitter;

    public:
        FluidSystem(double radius, double density, size_t size, bool fill = true);

        void updateNormals();        
        void updatePressures(double stiffness);
        void updateDensities(const std::vector<BoundarySystem> &boundaries);

        Eigen::Vector3d pressureAccFluid(const size_t i, const size_t j, const double ratio_i, const double ratio_j);
        Eigen::Vector3d viscAccFluid(const size_t i, const size_t j);
        Eigen::Vector3d pressureAccBoundary(const size_t i, const size_t k, const double ratio, const BoundarySystem& boundary);
        Eigen::Vector3d viscAccBoundary(const size_t i, const size_t k, const BoundarySystem& boundary);
        Eigen::Vector3d tensionForce(const size_t i, const size_t j);
        Eigen::Vector3d cohesionForce(const size_t i, const size_t j);
        Eigen::Vector3d curvatureForce(const size_t i, const size_t j);
        Eigen::Vector3d adhesionForce(const size_t i, const size_t j, const BoundarySystem& boundary);

        Eigen::Vector3d smoothingTerm(const size_t i, const size_t j);


        // Setter & Getter
        void addToParticleVel(const size_t i, const Eigen::Vector3d vel) { m_velocities[i] += vel; }
        void addToParticlePos(const size_t i, const Eigen::Vector3d pos) { m_positions[i] += pos; }
        void addToParticleAcc(const size_t i, const Eigen::Vector3d acc) { m_accelerations[i] += acc; }
        void setParticleAcc(const size_t i, const Eigen::Vector3d acc) { m_accelerations[i] = acc; }

        void setPrevPos(std::vector<Eigen::Vector3d> &prevPos) { m_prevPositions = prevPos; }
        std::vector<Eigen::Vector3d>& getPrevPos() { return m_prevPositions; }
        
        double getParticleDensity(size_t i) const { return m_densities[i]; }
        double getParticlePressure(size_t i) const { return m_pressures[i]; }
        Eigen::Vector3d getParticleNormal(size_t i) const { return m_normals[i]; }
        Eigen::Vector3d getParticleAcc(size_t i) const { return m_accelerations[i]; }
        Eigen::Vector3d getParticlePrevPos(size_t i) const { return m_prevPositions[i]; }

        const std::vector<double> &getDensities() const { return m_densities; }
        const std::vector<double> &getPressures() const { return m_pressures; }        
        const std::vector<Eigen::Vector3d>& getNormals() const { return m_normals; }
        const std::vector<Eigen::Vector3d>& getAccelerations() const { return m_accelerations; }                
        
        double getC() const { return m_c; }
        double getGamma() const { return m_gamma; }

        void setGamma(const double gamma) { m_gamma = gamma; }
      private:            
        Eigen::Vector3d normal(const size_t i);

        std::vector<double> m_pressures; // last updated particle pressures
        std::vector<double> m_densities; // last updated particle densities
        std::vector<Vector3d> m_normals; // normals for surface tension calculations
        std::vector<Eigen::Vector3d> m_accelerations; // acummulated particle accelerations
        std::vector<Eigen::Vector3d> m_prevPositions; // previous positions;
        
        double m_c = 0.3; // If particlesize is 0.1
        double m_gamma = 1.0;
        Kernel::Cohesion::Table m_cohesionWeightLookup;
        Kernel::Adhesion::Table m_adhesionWeightLookup;
    };
} // namespace learnSPH::System
