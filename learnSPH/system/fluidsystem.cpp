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
    m_densities.resize(size);
    m_pressures.resize(size);
    m_accelerations.resize(size);
    m_normalizedDensities.resize(size);
    m_normals.resize(size);
    m_tensionForces.resize(size);
    m_adhesionForces.resize(size);

    
    if (fill) {
        std::fill(m_densities.begin(), m_densities.end(), 0.0);
        std::fill(m_pressures.begin(), m_pressures.end(), 0.0);
        std::fill(m_accelerations.begin(), m_accelerations.end(), Vector3d(0.0, 0.0, 0.0));
        std::fill(m_normalizedDensities.begin(), m_normalizedDensities.end(), 0.0);
        std::fill(m_normals.begin(), m_normals.end(), Vector3d(0.0, 0.0, 0.0));
        std::fill(m_tensionForces.begin(), m_tensionForces.end(), Vector3d(0.0, 0.0, 0.0));
        std::fill(m_adhesionForces.begin(), m_adhesionForces.end(), Vector3d(0.0, 0.0, 0.0));
    }
}

void FluidSystem::clearTensionForces()
{
    std::fill(m_tensionForces.begin(), m_tensionForces.end(), Vector3d(0.0, 0.0, 0.0));
}

void FluidSystem::clearAdhesionForces()
{
    std::fill(m_adhesionForces.begin(), m_adhesionForces.end(), Vector3d(0.0, 0.0, 0.0));
}

void FluidSystem::addParticleTensionForce(size_t i, Vector3d force)
{
    m_tensionForces[i] += force;
}

void FluidSystem::addParticleAdhesionForce(size_t i, Vector3d force)
{
    m_adhesionForces[i] += force;
}

void FluidSystem::updatePressures(double stiffness)
{
    //#pragma omp parallel for
    for (size_t i = 0; i < getSize(); i++) {
        m_pressures[i] = std::max(0.0, stiffness * (m_densities[i] - m_restDensity));
    }
}

void FluidSystem::updateDensities(const std::vector<BoundarySystem> &boundaries)
{
    // get neighborhood information of fluid particle point set
    CompactNSearch::PointSet const& fluidPS = mp_nsearch->point_set(m_pointSetID);

    // iterate fluid particles
    //#pragma omp parallel for
    for (size_t i = 0; i < fluidPS.n_points(); i++) {
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


void FluidSystem::updateAccelerations(const std::vector<BoundarySystem> &boundaries,
                                      bool pressure, bool viscosity, bool external, bool tension, bool adhesion)
{
    //#pragma omp parallel for
    for (size_t i = 0; i < getSize(); i++) {
        m_accelerations[i] = Eigen::Vector3d(0.0, 0.0, 0.0);
    }
    if (pressure) {
      for (size_t i = 0; i < getSize(); i++) 
          m_accelerations[i] += particlePressureAcc(i, boundaries);      
    }
    if (viscosity) {
      for (size_t i = 0; i < getSize(); i++) 
          m_accelerations[i] += particleViscosityAcc(i, boundaries);
    }
    if (external) {
      for (size_t i = 0; i < getSize(); i++) 
          m_accelerations[i] += m_forces[i] / m_particleMass; // Check if this is actually correct, division should be by mass
    }
    if (tension) {
        for (size_t i = 0; i < getSize(); i++)
            m_accelerations[i] += m_tensionForces[i] / m_particleMass; // Check if this is actually correct, division should be by mass
    }

    if (adhesion) {
        for (size_t i = 0; i < getSize(); i++)
            m_accelerations[i] += m_adhesionForces[i] / m_particleMass; // Check if this is actually correct, division should be by mass
    }
}


Vector3d FluidSystem::particlePressureAcc(size_t i, const std::vector<BoundarySystem> &boundaries)
{
    const size_t id = getPointSetID();
    const Eigen::Vector3d &pos = getParticlePos(i);
    const double pressureDensityRatio = m_pressures[i] / pow(m_densities[i], 2);
   
    CompactNSearch::PointSet const& fluidPS = mp_nsearch->point_set(id);
    
    // Fluid contribution
    Eigen::Vector3d fluidContrib(0.0, 0.0, 0.0);
    fluidContrib += m_particleMass * 2 * pressureDensityRatio
        * m_kernelLookup.gradWeight(pos, pos);
    for (size_t j = 0; j < fluidPS.n_neighbors(id, i); j++) {
        const unsigned int k = fluidPS.neighbor(id, i, j);

        double pressureTerm = pressureDensityRatio
            + (m_pressures[k] / pow(m_densities[k], 2));
        
        fluidContrib += m_particleMass * pressureTerm
            * m_kernelLookup.gradWeight(pos, m_positions[k]);
    }

    // Boundary contribution
    Eigen::Vector3d boundaryContrib(0.0, 0.0, 0.0);
    for (const BoundarySystem &boundary : boundaries) {
        const size_t boundaryID = boundary.getPointSetID();
        for (size_t j = 0; j < fluidPS.n_neighbors(boundaryID, i); j++) {
            const unsigned int k = fluidPS.neighbor(boundaryID, i, j);
            
            boundaryContrib += m_restDensity * boundary.getParticleVolume(k)
                * pressureDensityRatio
                * m_kernelLookup.gradWeight(pos, boundary.getParticlePos(k));
        }
    }
    
    return -(fluidContrib + boundaryContrib);
}

Vector3d FluidSystem::particleViscosityAcc(size_t i, const std::vector<BoundarySystem> &boundaries)
{
    const size_t id = getPointSetID();
    const Eigen::Vector3d &pos = getParticlePos(i);
   
    CompactNSearch::PointSet const& fluidPS = mp_nsearch->point_set(id);
    
    // Fluid contribution
    Eigen::Vector3d fluidContrib (0.0, 0.0, 0.0);
    if (getViscosity() >= 0.0) {
        for (size_t j = 0; j < fluidPS.n_neighbors(id, i); j++) {
            const unsigned int k = fluidPS.neighbor(id, i, j);

            Eigen::Vector3d posDiff = getParticlePos(i) - getParticlePos(k);
            fluidContrib +=
                (m_particleMass / getParticleDensity(k)) *
                (getParticleVel(i) - getParticleVel(k)) *
                posDiff.dot(m_kernelLookup.gradWeight(pos, getParticlePos(k))) /
                (pow(posDiff.norm(), 2) + 0.01 * pow(m_smoothingLength, 2));
        }
        fluidContrib *= getViscosity();
    }

    // Boundary contribution
    Eigen::Vector3d boundaryContrib(0.0, 0.0, 0.0);
    for (const BoundarySystem& boundary : boundaries) {
        if (boundary.getViscosity() == 0.0)
            continue;
        
        const size_t boundaryID = boundary.getPointSetID();
        for (size_t j = 0; j < fluidPS.n_neighbors(boundaryID, i); j++) {
            const unsigned int k = fluidPS.neighbor(boundaryID, i, j);

            Eigen::Vector3d posDiff = pos - boundary.getParticlePos(k);
            
            boundaryContrib += boundary.getViscosity() 
                * boundary.getParticleVolume(k) * getParticleVel(i)
                * posDiff.dot(m_kernelLookup.gradWeight(pos, boundary.getParticlePos(k)))
                * (pow(posDiff.norm(), 2) + 0.01 * pow(m_smoothingLength, 2)); // should be / in the beginning of the line
        }
    }
    
    return 2.0 * (fluidContrib + boundaryContrib);
}

void FluidSystem::updateNormals() {
    for (size_t i = 0; i < getSize(); i++) {
        m_normals[i] = normal(i);
    }
}

Eigen::Vector3d FluidSystem::normal(const size_t i) {
    Eigen::Vector3d normal = (m_particleMass / getParticleDensity(i)) * m_kernelLookup
        .gradWeight(getParticlePos(i), getParticlePos(i));
    CompactNSearch::PointSet const& fluidPS = mp_nsearch->point_set(m_pointSetID);
    for (size_t j = 0; j < fluidPS.n_neighbors(m_pointSetID, i); j++) {
        const unsigned int pid = fluidPS.neighbor(m_pointSetID, i, j);
        normal += (m_particleMass / getParticleDensity(pid)) * m_kernelLookup
            .gradWeight(getParticlePos(i), getParticlePos(pid));
    }

    return m_c * normal;
}

Eigen::Vector3d FluidSystem::pressureAccFluid(const size_t i, const size_t j, const double ratio_i, const double ratio_j) {
    return m_particleMass * (ratio_i + ratio_j) * m_kernelLookup.gradWeight(m_positions[i], m_positions[j]);
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
    return m_restDensity * boundary.getParticleVolume(k) * ratio * m_kernelLookup.gradWeight(m_positions[i], m_positions[k]);
}

Eigen::Vector3d FluidSystem::viscAccBoundary(const size_t i, const size_t k, const BoundarySystem& boundary) {
    const Eigen::Vector3d posDiff = m_positions[i] - m_positions[k];
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

    return -m_gamma * m_particleMass * m_particleMass * cohesionWeight(diffNorm) * diff.normalized();
}

Eigen::Vector3d FluidSystem::curvatureForce(const size_t i, const size_t j) {
    return -m_gamma * m_particleMass * (m_normals[i] - m_normals[j]);
}

Eigen::Vector3d FluidSystem::adhesionForce(const size_t i, const size_t k, const BoundarySystem& boundary) {
    const Eigen::Vector3d diff = m_positions[i] - boundary.getParticlePos(k);
    const double diffNorm = diff.norm();
    return -boundary.getBeta() * m_particleMass * m_particleMass * adhesionWeight(diffNorm)* diff.normalized();
}

Eigen::Vector3d FluidSystem::smoothingTerm(const size_t i, const size_t j) {
    return 2.0 * m_particleMass * (m_velocities[j] - m_velocities[i]) / (m_densities[i] + m_densities[j]) 
        * m_kernelLookup.weight(m_positions[i], m_positions[j]);
}

double FluidSystem::cohesionWeight(const double r) {
    const double alpha = 32.0 / (M_PI * m_c * m_c * m_c * m_c * m_c * m_c * m_c * m_c * m_c);
    if (0.0 <= r && r <= m_c / 2.0) {
        return alpha * 2 * (m_c - r) * (m_c - r) * (m_c - r) * r * r * r - (m_c * m_c * m_c * m_c * m_c * m_c) / 64.0;
    }
    else if (m_c / 2.0 < r && r <= m_c) {
        return alpha * (m_c - r) * (m_c - r) * (m_c - r) * r * r * r;
    }
    return 0.0;
}

double FluidSystem::adhesionWeight(const double r) {
    const double alpha = 0.007 / pow(m_c, 3.25);

    if (m_c / 2.0 <= r && r <= m_c) {
        return alpha * pow((-4.0 * r * r) / m_c + 6.0 * r - 2.0 * m_c, (1.0 / 4.0));
    }
    return 0;
}
