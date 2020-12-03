#include "fluidsystem.h"
#include "boundarysystem.h"
#include <algorithm>

using namespace learnSPH;
using namespace learnSPH::System;
using namespace learnSPH::Kernel;
using namespace CompactNSearch;

FluidSystem::FluidSystem(double radius, double density, size_t size, bool fill)
    : ParticleSystem(radius, density, size, fill)
{
    m_kernelLookup.generateTable(m_smoothingLength, 1000);
    
    m_densities.resize(size);
    m_pressures.resize(size);
    m_accelerations.resize(size);
    
    if (fill) {
        std::fill(m_densities.begin(), m_densities.end(), 0.0);
        std::fill(m_pressures.begin(), m_pressures.end(), 0.0);
        std::fill(m_accelerations.begin(), m_accelerations.end(), Vector3d(0.0, 0.0, 0.0));
    }
}

void FluidSystem::updatePressures(double stiffness)
{
    for (size_t i = 0; i < getSize(); i++) {
        m_pressures[i] = std::max(0.0, stiffness * (m_densities[i] - m_restDensity));
    }
}

void FluidSystem::updateDensities(const std::vector<BoundarySystem> &boundaries)
{
    // get neighborhood information of fluid particle point set
    CompactNSearch::PointSet const& fluidPS = mp_nsearch->point_set(m_pointSetID);
    m_densities.resize(fluidPS.n_points());
    
    // iterate fluid particles
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
            size_t n_neighbors = fluidPS.n_neighbors(boundaryID, i);
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


void FluidSystem::updateAccelerations(const std::vector<BoundarySystem> &boundaries)
{
    for (size_t i = 0; i < getSize(); i++) {
        Eigen::Vector3d pressureAcc = particlePressureAcc(i, boundaries);
        Eigen::Vector3d viscosityAcc = particleViscosityAcc(i, boundaries);
        Eigen::Vector3d externalAcc = m_forces[i] / m_restDensity;

        m_accelerations[i] = pressureAcc + viscosityAcc + externalAcc;
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
                * (pow(posDiff.norm(), 2) + 0.01 * pow(m_smoothingLength, 2));
        }
    }
    

    return 2.0 * (fluidContrib + boundaryContrib);
}

