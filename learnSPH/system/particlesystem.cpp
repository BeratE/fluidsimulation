#include "particlesystem.h"

using namespace learnSPH;
using namespace learnSPH::System;
using namespace CompactNSearch;

ParticleSystem::ParticleSystem(size_t size, bool fill)
{
    m_positions.resize(size);
    m_velocities.resize(size);
    m_accelerations.resize(size);
    m_forces.resize(size);
    if (fill) {
        std::fill(m_positions.begin(), m_positions.end(),
                  Eigen::Vector3d(0.0, 0.0, 0.0));
        std::fill(m_velocities.begin(), m_velocities.end(),
                  Eigen::Vector3d(0.0, 0.0, 0.0));
        std::fill(m_accelerations.begin(), m_accelerations.end(),
                  Eigen::Vector3d(0.0, 0.0, 0.0));
        std::fill(m_forces.begin(), m_forces.end(),
                  Eigen::Vector3d(0.0, 0.0, 0.0));
    }
}

double ParticleSystem::smoothingLength() const
{
    return 2 * m_particleRadius * Kernel::Parameter::TUNING;
}

double ParticleSystem::particleMass() const
{
    return pow(2 * m_particleRadius, 3) * m_restDensity;
}

void ParticleSystem::addToNeighborhood(const std::shared_ptr<NeighborhoodSearch> &nsearch)
{
    mp_nsearch = nsearch;
    m_pointSetID = mp_nsearch->add_point_set(m_positions.front().data(),
                                             m_positions.size());
}

void ParticleSystem::addParticlePos(size_t i, Eigen::Vector3d pos)
{
    m_positions[i] += pos;
}

void ParticleSystem::addParticleVel(size_t i, Eigen::Vector3d vel)
{
    m_velocities[i] += vel;
}

void ParticleSystem::addParticleAcc(size_t i, Eigen::Vector3d acc)
{
    m_accelerations[i] += acc;
}

void ParticleSystem::addParticleForce(size_t i, Eigen::Vector3d force)
{
    m_forces[i] += force;
}

void ParticleSystem::clearForces()
{
    std::fill(m_forces.begin(), m_forces.end(),
              Eigen::Vector3d(0.0, 0.0, 0.0));
}