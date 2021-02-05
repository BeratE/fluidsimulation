#pragma once
#include "../kernel.h"
#include <vector>
#include <memory>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <CompactNSearch/CompactNSearch.h>

using namespace Eigen;
using namespace CompactNSearch;

namespace learnSPH::System {
    class ParticleSystem {
        friend class ParticleEmitter;

    public:
        ParticleSystem(double radius, double density,
                       size_t size, bool fill = true);

        void addToNeighborhood(const std::shared_ptr<NeighborhoodSearch> &nsearch);

        void addParticlePos(size_t i, Vector3d pos);
        void addParticleVel(size_t i, Vector3d vel);

        // Setter & Getter
        void setViscosity(double value) { m_viscosity = value; }
        void setPointSetID(double id) { m_pointSetID = id; }

        void setParticlePos(size_t i, Vector3d pos) {m_positions[i] = pos;}
        void setParticleVel(size_t i, Vector3d vel) {m_velocities[i] = vel;}        
        
        size_t getSize() const { return m_positions.size(); }
        size_t getPointSetID() const { return m_pointSetID; }
        double getViscosity() const {return m_viscosity; }
        double getRestDensity() const { return m_restDensity; }
        double getParticleMass() const { return m_particleMass; }
        double getParticleRadius() const { return m_particleRadius; }
        double getSmoothingLength() const { return m_smoothingLength; }

        Eigen::Vector3d getParticlePos(size_t i) const {return m_positions[i];}
        Eigen::Vector3d getParticleVel(size_t i) const {return m_velocities[i];}

        const std::vector<Vector3d>& getPositions() const { return m_positions; }
        const std::vector<Vector3d>& getVelocities() const { return m_velocities; }

        const Kernel::CubicSpline::Table getKernelLookUp() const { return m_kernelLookup;  }
        const std::shared_ptr<NeighborhoodSearch> getNSearch() const { return mp_nsearch;  }

    protected:
        long m_pointSetID = -1; // CompactNSearch PointSet ID
        double m_viscosity = 0.0; // Viscosity, different for fluid and boundary
        double m_restDensity = 1000.0;
        double m_particleMass = 1.0;
        double m_particleRadius = 1.0; // (cubic root of volume) / 2
        double m_smoothingLength = 1.0;
        
        std::vector<Vector3d> m_positions; // Particle Positions
        std::vector<Vector3d> m_velocities; // Particle Velocities        

        std::shared_ptr<NeighborhoodSearch> mp_nsearch; // reference to neighborhood information

        Kernel::CubicSpline::Table m_kernelLookup;
    };
} // namespace learnSPH::System
