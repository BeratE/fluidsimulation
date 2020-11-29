#pragma once
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <CompactNSearch/CompactNSearch.h>

namespace learnSPH {
    class ParticleSystem {
        friend class Emitter;

    protected: ParticleSystem() {}
        
    public:
        static constexpr double B = 5000.0;

        double smoothingLength();
        size_t getSize() { return m_positions.size(); }
        double getRestDensity() { return m_restDensity; }
        double getParticleMass() { return m_particleMass; }
        double getParticleRadius() { return m_particleRadius; }
        double getPointSetID() { return m_pointSetID; }
        Eigen::Vector3d getParticlePosition(size_t i) {return m_positions[i];}
        std::vector<Eigen::Vector3d> getPositions() { return m_positions; }

        void setPointSetID(double id) { m_pointSetID = id; }
        void setParticlePosition(size_t i, Eigen::Vector3d pos) {m_positions[i] = pos;}

      protected:
        long m_pointSetID = -1; // CompactNSearch PointSet ID
        double m_restDensity = 1.0;
        double m_particleMass = 1.0;   // Constant particle density for all particles
        double m_particleRadius = 1.0; // (cubic root of volume) / 2
        double m_viscosity = 0.0; // Viscosity, different for fluid and boundary
        std::vector<Eigen::Vector3d> m_positions; // Particle Positions
    };
    
} // namespace learnSPH



// /* Sample positions in a triangle specified by 3d points a, b , c with
//  * sampling distance apart from each other (hexagonal).
//  * @param a, b , c         - 3d points of the triangle.
//  * @param samplingDistance - distance the sample points are apart.
//  * @return Vector of sampled 3d points within the triangle. */
// std::vector<Eigen::Vector3d> samplePositionsTriangle(Eigen::Vector3d a,
//                                                      Eigen::Vector3d b,
//                                                      Eigen::Vector3d c,
//                                                      double samplingDistance);

// /* Sample a hollow box spanning from the bottom left to the top right corner.
//  * @param bottomLeft - bottom left box corner.
//  * @param topRight   - top right box corner diagonal to bottom left.
//  * @param samplingDistance - particles are sampling distance apart.
//  * @return vector of positions of sample points. */
// std::vector<Eigen::Vector3d> samplePositionsBox(Eigen::Vector3d bottomLeft,
//                                                 Eigen::Vector3d topRight,
//                                                 double samplingDistance);

// /* Construct a boundary from the given positions, rest density and mass */
// BoundarySystem createBoundary(const std::vector<Eigen::Vector3d> &positions,
//                               double restDensity, double particleMass,
//                               double particleRadius);

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
