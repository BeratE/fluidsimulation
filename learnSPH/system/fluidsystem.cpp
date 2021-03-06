#define _USE_MATH_DEFINES
#include "fluidsystem.h"
#include "boundarysystem.h"
#include <math.h>
#include <algorithm>

using namespace learnSPH;
using namespace learnSPH::System;
using namespace learnSPH::Kernel;
using namespace CompactNSearch;

FluidSystem::FluidSystem(double radius, double density, size_t size, bool fill)
    : ParticleSystem(radius, density, size, fill)
{   
    m_c = 2 * m_smoothingLength;
    m_cohesionWeightLookup.generateTable(m_c, 1000);
    m_adhesionWeightLookup.generateTable(m_c, 1000);

    m_normals.resize(size);
    m_densities.resize(size);
    m_pressures.resize(size);
    m_accelerations.resize(size);
    
    m_prevPositions = m_positions;
    
    if (fill) {
        std::fill(m_densities.begin(), m_densities.end(), 0.0);
        std::fill(m_pressures.begin(), m_pressures.end(), 0.0);
        std::fill(m_normals.begin(), m_normals.end(), Eigen::Vector3d(0.0, 0.0, 0.0));
        std::fill(m_accelerations.begin(), m_accelerations.end(), Eigen::Vector3d(0.0, 0.0, 0.0));
    }
}

void FluidSystem::updatePressures(double stiffness)
{
    #pragma omp for schedule(static)
    for (int i = 0; i < getSize(); i++) {
        m_pressures[i] = std::max(0.0, stiffness * (m_densities[i] - m_restDensity));
    }
}

void FluidSystem::updateDensities(const std::vector<BoundarySystem> &boundaries)
{
    // get neighborhood information of fluid particle point set
    CompactNSearch::PointSet const& fluidPS = mp_nsearch->point_set(m_pointSetID);

    // iterate fluid particles
    #pragma omp for schedule(static) 
    for (int i = 0; i < fluidPS.n_points(); i++) {
        const Eigen::Vector3d &fpPos = m_positions[i];

        // Fluid contribution
        double fluidDensity = 0.0;
        fluidDensity += m_kernelLookup.weight(fpPos, fpPos);        
        for (size_t j = 0; j < fluidPS.n_neighbors(m_pointSetID, i); j++) {
            const unsigned int k = fluidPS.neighbor(m_pointSetID, i, j);
            fluidDensity += m_kernelLookup.weight(fpPos, m_positions[k]);
        }
        fluidDensity *= m_particleMass;

        // Boundary contribution
        double boundaryDensity = 0.0;
        for (const BoundarySystem &boundary : boundaries) {
            double density = 0.0;
            const size_t boundaryID = boundary.getPointSetID();
            const size_t n_neighbors = fluidPS.n_neighbors(boundaryID, i);
            for (size_t j = 0; j < n_neighbors; j++) {
                const unsigned int k = fluidPS.neighbor(boundaryID, i, j);
                density += boundary.getParticleVolume(k)
                    * m_kernelLookup.weight(fpPos, boundary.getParticlePos(k));
            }
            density *= boundary.getRestDensity();
            boundaryDensity += density;
        }

        m_densities[i] = fluidDensity + boundaryDensity;
    }
}

void FluidSystem::updateNormals() {
    #pragma omp for schedule(static) 
    for (int i = 0; i < getSize(); i++) {
        m_normals[i] = normal(i);
    }
}

Eigen::Vector3d FluidSystem::normal(const size_t i) {
    Eigen::Vector3d normal = (m_particleMass / getParticleDensity(i))
        * m_kernelLookup.gradWeight(getParticlePos(i), getParticlePos(i));
    
    CompactNSearch::PointSet const& fluidPS = mp_nsearch->point_set(m_pointSetID);
    for (size_t j = 0; j < fluidPS.n_neighbors(m_pointSetID, i); j++) {
        const unsigned int pid = fluidPS.neighbor(m_pointSetID, i, j);
        normal += (m_particleMass / getParticleDensity(pid))
            * m_kernelLookup.gradWeight(getParticlePos(i), getParticlePos(pid));
    }

    return m_c * normal;
}



Eigen::Vector3d FluidSystem::pressureAccFluid(const size_t i, const size_t j, const double ratio_i, const double ratio_j) {
    return m_particleMass * (ratio_i + ratio_j) *
        m_kernelLookup.gradWeight(m_positions[i], m_positions[j]);
}

Eigen::Vector3d FluidSystem::viscAccFluid(const size_t i, const size_t j) {
    const Eigen::Vector3d posDiff = m_positions[i] - m_positions[j];
    const double posDiffNorm = posDiff.norm();
    
    const Eigen::Vector3d velDiff = m_velocities[i] - m_velocities[j];

    return 2.0 * m_viscosity * (m_particleMass / m_densities[j]) * velDiff 
        * (posDiff.dot(m_kernelLookup.gradWeight(m_positions[i], m_positions[j]))) 
        / (posDiffNorm * posDiffNorm + 0.01 * m_smoothingLength * m_smoothingLength);
}

Eigen::Vector3d FluidSystem::pressureAccBoundary(const size_t i, const size_t k, const double ratio, const BoundarySystem& boundary) {
    return m_restDensity * boundary.getParticleVolume(k) * ratio *
        m_kernelLookup.gradWeight(m_positions[i], boundary.getParticlePos(k));
}

Eigen::Vector3d FluidSystem::viscAccBoundary(const size_t i, const size_t k, const BoundarySystem& boundary) {
    const Eigen::Vector3d posDiff = m_positions[i] - boundary.getParticlePos(k);
    const double posDiffNorm = posDiff.norm();
    return 2.0 * boundary.getViscosity()
        * boundary.getParticleVolume(k) * getParticleVel(i)
        * posDiff.dot(m_kernelLookup.gradWeight(m_positions[i], boundary.getParticlePos(k)))
        / (posDiffNorm * posDiffNorm + 0.01 * m_smoothingLength * m_smoothingLength);
}

Eigen::Vector3d FluidSystem::tensionForce(const size_t i, const size_t j) {
    const double K = (2.0 * m_restDensity) / (m_densities[i] + m_densities[j]);
    return K * (cohesionForce(i, j) + curvatureForce(i, j));
}

Eigen::Vector3d FluidSystem::cohesionForce(const size_t i, const size_t j) {
    const Eigen::Vector3d diff = m_positions[i] - m_positions[j];
    const double diffNorm = diff.norm();
    return -m_gamma * m_particleMass * m_particleMass
        * Cohesion::weight(diffNorm, m_c) * diff.normalized();
}

Eigen::Vector3d FluidSystem::curvatureForce(const size_t i, const size_t j) {
    return -m_gamma * m_particleMass * (m_normals[i] - m_normals[j]);
}

Eigen::Vector3d FluidSystem::adhesionForce(const size_t i, const size_t k,
                                           const BoundarySystem& boundary) {
    const Eigen::Vector3d diff = m_positions[i] - boundary.getParticlePos(k);
    const double diffNorm = diff.norm();
    return -boundary.getBeta() * m_particleMass * boundary.getParticleMass(k)
        * Adhesion::weight(diffNorm, m_c)* diff.normalized();
}

Eigen::Vector3d FluidSystem::smoothingTerm(const size_t i, const size_t j) {
    return 2.0 * m_particleMass * (m_velocities[j] - m_velocities[i]) / (m_densities[i] + m_densities[j]) 
        * m_kernelLookup.weight(m_positions[i], m_positions[j]);
}
