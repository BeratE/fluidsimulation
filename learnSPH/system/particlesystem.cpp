#include "particlesystem.h"
#include <float.h>

using namespace learnSPH;
using namespace learnSPH::System;

ParticleSystem::ParticleSystem(double radius, double density, size_t size, bool fill)
    : m_particleRadius(radius), m_restDensity(density)
{
    m_smoothingLength = 2 * m_particleRadius * Kernel::Parameter::TUNING;
    std::cout << 2 * m_smoothingLength;
    m_kernelLookup.generateTable(m_smoothingLength, 1000);
    m_particleMass = pow(2 * m_particleRadius, 3) * m_restDensity;
    
    m_positions.resize(size);
    m_velocities.resize(size);
    m_normalizedDensities.resize(size);
    
    if (fill) {
        std::fill(m_positions.begin(), m_positions.end(), Vector3d(0.0, 0.0, 0.0));
        std::fill(m_velocities.begin(), m_velocities.end(), Vector3d(0.0, 0.0, 0.0));
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



void ParticleSystem::updateNormalizedDensities() {
    CompactNSearch::PointSet const ps = mp_nsearch->point_set(m_pointSetID);
    for (int fpIdx = 0; fpIdx < ps.n_points(); fpIdx++) {
        double normalizedDensity = m_kernelLookup.weight(m_positions[fpIdx], m_positions[fpIdx]);
        for (size_t fnId = 0; fnId < ps.n_neighbors(m_pointSetID, fpIdx); fnId++) {
            const unsigned int nId = ps.neighbor(m_pointSetID, fpIdx, fnId);
            normalizedDensity += m_kernelLookup.weight(m_positions[fpIdx], m_positions[nId]);
        }
        m_normalizedDensities[fpIdx] = normalizedDensity;
    }
}
