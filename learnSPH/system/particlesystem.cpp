#include "particlesystem.h"
#include <float.h>

using namespace learnSPH;
using namespace learnSPH::System;

ParticleSystem::ParticleSystem(double radius, double density, size_t size, bool fill)
    : m_particleRadius(radius), m_restDensity(density)
{
    m_smoothingLength = 2 * m_particleRadius * Kernel::Parameter::TUNING;
    m_particleMass = pow(2 * m_particleRadius, 3) * m_restDensity;
    
    m_positions.resize(size);
    m_velocities.resize(size);
    m_forces.resize(size);
    
    if (fill) {
        std::fill(m_positions.begin(), m_positions.end(), Vector3d(0.0, 0.0, 0.0));
        std::fill(m_velocities.begin(), m_velocities.end(), Vector3d(0.0, 0.0, 0.0));
        std::fill(m_forces.begin(), m_forces.end(), Vector3d(0.0, 0.0, 0.0));
    }
}

void ParticleSystem::addToNeighborhood(const std::shared_ptr<NeighborhoodSearch> &nsearch)
{
    mp_nsearch = nsearch;
    m_pointSetID = mp_nsearch->add_point_set(m_positions.front().data(),
                                             m_positions.size());
}

void ParticleSystem::addParticlePos(size_t i, Vector3d pos)
{
    m_positions[i] += pos;
}

void ParticleSystem::addParticleVel(size_t i, Vector3d vel)
{
    m_velocities[i] += vel;
}


void ParticleSystem::addParticleForce(size_t i, Vector3d force)
{
    m_forces[i] += force;
}

void ParticleSystem::clearForces()
{
    std::fill(m_forces.begin(), m_forces.end(), Vector3d(0.0, 0.0, 0.0));
}

void ParticleSystem::boundingBox(Eigen::Vector3d& bottomLeft, Eigen::Vector3d& upperRight) const {
    double minX = DBL_MAX, minY = DBL_MAX, minZ = DBL_MAX;
    double maxX = DBL_MIN, maxY = DBL_MIN, maxZ = DBL_MIN;
    for (Eigen::Vector3d position : m_positions) {
        if (position.x() < minX)
            minX = position.x();
        if (position.y() < minY)
            minY = position.y();
        if (position.z() < minZ)
            minZ = position.z();
        if (position.x() > maxX)
            maxX = position.x();
        if (position.y() > maxY)
            maxY = position.y();
        if (position.z() > maxZ)
            maxZ = position.z();
    }
    bottomLeft(0) = minX;
    bottomLeft(1) = minY;
    bottomLeft(2) = minZ;
    upperRight(0) = maxX;
    upperRight(1) = maxY;
    upperRight(2) = maxZ;
    
}
