#pragma once
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <CompactNSearch/CompactNSearch.h>


namespace learnSPH {
    class ParticleSystem {
        friend class ParticleEmitter;

    protected:
        ParticleSystem() {}
        
    public:
        static constexpr double B = 5000.0;        

        double smoothingLength() const;
        double particleMass() const;
        
        void addToNeighborhood(CompactNSearch::NeighborhoodSearch &nsearch);

        // Setter & Getter
        size_t getSize() const { return m_positions.size(); }
        double getRestDensity() const { return m_restDensity; }
        double getParticleRadius() const { return m_particleRadius; }
       
        long getPointSetID() const { return m_pointSetID; }
        void setPointSetID(double id) { m_pointSetID = id; }

        Eigen::Vector3d getParticlePos(size_t i) const {return m_positions[i];}
        void setParticlePos(size_t i, Eigen::Vector3d pos) {m_positions[i] = pos;}

        const std::vector<Eigen::Vector3d>& getPositions() const { return m_positions; }

    protected:
        long m_pointSetID = -1; // CompactNSearch PointSet ID
        double m_restDensity = 1.0;
        double m_particleRadius = 1.0; // (cubic root of volume) / 2
        double m_viscosity = 0.0; // Viscosity, different for fluid and boundary
        std::vector<Eigen::Vector3d> m_positions; // Particle Positions
    };
    
} // namespace learnSPH









// /* Estimates the positions of particles at desiredTime by interpolating the
//  * positions at prevTime and curTime
//  * @param &prevPositions - Previous positions of the particles
//  * @param &curPositions - Current positions of the particles
//  * @param prevTime - Time at which the prevPositions were determined
//  * @param curTime - Time at which the curPositions were determined
//  * @param desiredTime - Time at which the positions of the particles should be
//  * estimated. Requires prevTime <= desiredTime <= curTime.
//  * @returns interpolated positions of particles at desiredTime */
// std::vector<Eigen::Vector3d>
// interpolatePositions(const std::vector<Eigen::Vector3d> &prevPositions,
//                      const std::vector<Eigen::Vector3d> &curPositions,
//                      const double prevTime, const double curTime,
//                      const double desiredTime);

// double calculatePressure(const double particleDensity,
//                          const double restDensity);
