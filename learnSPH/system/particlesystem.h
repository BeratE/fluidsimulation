#pragma once
#include "../kernel.h"
#include <vector>
#include <memory>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <CompactNSearch/CompactNSearch.h>

namespace learnSPH::System {
    class ParticleSystem {
        friend class ParticleEmitter;

    public:
        ParticleSystem() {}
        ParticleSystem(size_t size, bool fill = true);

        double smoothingLength() const;
        double particleMass() const;
        
        void addToNeighborhood(const std::shared_ptr<CompactNSearch::NeighborhoodSearch> &nsearch);
      
        void clearForces();
        void addParticlePos(size_t i, Eigen::Vector3d pos);
        void addParticleVel(size_t i, Eigen::Vector3d vel);
        void addParticleAcc(size_t i, Eigen::Vector3d acc);
        void addParticleForce(size_t i, Eigen::Vector3d force);
        
        // Setter & Getter
        size_t getSize() const { return m_positions.size(); }
        double getViscosity() const {return m_viscosity; }
        double getRestDensity() const { return m_restDensity; }
        double getParticleRadius() const { return m_particleRadius; }
        size_t getPointSetID() const { return m_pointSetID; }

        void setViscosity(double value) {m_viscosity = value;}
        void setPointSetID(double id) { m_pointSetID = id; }

        Eigen::Vector3d getParticlePos(size_t i) const {return m_positions[i];}
        Eigen::Vector3d getParticleVel(size_t i) const {return m_velocities[i];}
        Eigen::Vector3d getParticleAcc(size_t i) const {return m_accelerations[i];}
        Eigen::Vector3d getParticleForce(size_t i) const {return m_forces[i];}
        void setParticlePos(size_t i, Eigen::Vector3d pos) {m_positions[i] = pos;}
        void setParticleVel(size_t i, Eigen::Vector3d vel) {m_velocities[i] = vel;}        
        void setParticleAcc(size_t i, Eigen::Vector3d acc) {m_accelerations[i] = acc;}

        const std::vector<Eigen::Vector3d>& getForces() const { return m_forces; }
        const std::vector<Eigen::Vector3d>& getPositions() const { return m_positions; }
        const std::vector<Eigen::Vector3d>& getVelocities() const { return m_velocities; }
        const std::vector<Eigen::Vector3d>& getAccelerations() const { return m_accelerations; }

    protected:
        long m_pointSetID = -1; // CompactNSearch PointSet ID
        double m_viscosity = 0.0; // Viscosity, different for fluid and boundary
        double m_restDensity = 1000.0;
        double m_particleRadius = 1.0; // (cubic root of volume) / 2
        std::vector<Eigen::Vector3d> m_positions; // Particle Positions
        std::vector<Eigen::Vector3d> m_velocities; // Particle Velocities
        std::vector<Eigen::Vector3d> m_accelerations; // Particle Accelerations
        std::vector<Eigen::Vector3d> m_forces; // Particle Force Accumulator

        std::shared_ptr<CompactNSearch::NeighborhoodSearch> mp_nsearch;
    };
} // namespace learnSPH::System
