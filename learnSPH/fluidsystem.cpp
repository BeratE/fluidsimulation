#include "fluidsystem.h"
#include "boundarysystem.h"
#include "kernel.h"
#include <algorithm>

using namespace learnSPH;
using namespace learnSPH::Kernel;
using namespace CompactNSearch;

FluidSystem::FluidSystem(size_t size, bool fill)
    : ParticleSystem(size, fill)
{
    m_densities.resize(size);
    m_pressures.resize(size);
    if (fill) {
        std::fill(m_densities.begin(), m_densities.end(), 0.0);
        std::fill(m_pressures.begin(), m_pressures.end(), 0.0);
    }
}

void FluidSystem::initTable()
{
    m_kernelTable.generateTable(smoothingLength(), 1000);
}

void FluidSystem::estimateDensity(NeighborhoodSearch &nsearch)
{
    estimateDensity(nsearch, std::vector<BoundarySystem>());
}

void FluidSystem::estimateDensity(NeighborhoodSearch &nsearch,
                                  const std::vector<BoundarySystem> &boundaries)
{
    // get neighborhood information of fluid particle point set
    CompactNSearch::PointSet const& fluidPS = nsearch.point_set(m_pointSetID);
    m_densities.resize(fluidPS.n_points());
    
    // iterate fluid particles
    for (size_t fpI = 0; fpI < fluidPS.n_points(); fpI++) {
        const Eigen::Vector3d &fpPos = m_positions[fpI];

        // Fluid contribution
        /* The current fluid particle itself is part of its own neighborhood. */
        double fluidDensity = 0.0;
        fluidDensity += m_kernelTable.weight(fpPos, fpPos);
        for (size_t fpN = 0; fpN < fluidPS.n_neighbors(m_pointSetID, fpI); fpN++) {
            const unsigned int fnI = fluidPS.neighbor(m_pointSetID, fpI, fpN);
            fluidDensity += CubicSpline::weight(fpPos, m_positions[fnI], smoothingLength());
        }
        fluidDensity *= particleMass();

        // Boundary contribution
        double boundaryDensity = 0.0;
        for (const BoundarySystem &boundary : boundaries) {
            double density = 0.0;
            size_t n_neighbors = fluidPS.n_neighbors(boundary.getPointSetID(), fpI);
            for (size_t bpI = 0; bpI < n_neighbors; bpI++) {
                const unsigned int bnI = fluidPS.neighbor(boundary.getPointSetID(),
                                                          fpI, bpI);
                density += boundary.getParticleVolume(bnI)
                    * CubicSpline::weight(fpPos, boundary.getParticlePos(bnI), smoothingLength());
            }
            density *= boundary.getRestDensity();
            boundaryDensity += density;
        }

        m_densities[fpI] = fluidDensity + boundaryDensity;
    }  
}

void FluidSystem::updatePressure()
{
    for (size_t i = 0; i < getSize(); i++) {
        m_pressures[i] = std::max(0.0, m_stiffness * (m_densities[i] - m_restDensity));
    }
}


void FluidSystem::updateAcceleration(CompactNSearch::NeighborhoodSearch& nsearch,
                                      const std::vector<BoundarySystem> &boundaries)
{
    for (size_t i = 0; i < getSize(); i++) {
        Eigen::Vector3d pressureAcc = particlePressureAcc(i, nsearch, boundaries);
        Eigen::Vector3d viscosityAcc = particleViscosityAcc(i, nsearch, boundaries);
        //Eigen::Vector3d externalAcc = m_forces[i] / m_restDensity;

        m_accelerations[i] = pressureAcc + viscosityAcc;
    }
}


Eigen::Vector3d FluidSystem::particlePressureAcc(size_t index,
                                                 CompactNSearch::NeighborhoodSearch& nsearch,
                                                 const std::vector<BoundarySystem> &boundaries)
{
    const size_t id = getPointSetID();
    const Eigen::Vector3d &pos = getParticlePos(index);
    const double particleMass = this->particleMass();
    const double particleDensityRatio = m_pressures[index] / pow(m_densities[index], 2);
   
    CompactNSearch::PointSet const& fluidPS = nsearch.point_set(id);
    
    // Fluid contribution
    Eigen::Vector3d fluidContrib(0.0, 0.0, 0.0);
    fluidContrib += particleMass * 2 * particleDensityRatio
        * CubicSpline::gradWeight(pos, pos, smoothingLength());
    for (size_t j = 0; j < fluidPS.n_neighbors(id, index); j++) {
        const unsigned int k = fluidPS.neighbor(id, index, j);

        double pressureTerm = particleDensityRatio
            + (m_pressures[k] / pow(m_densities[k], 2));
        
        fluidContrib += particleMass * pressureTerm
            * CubicSpline::gradWeight(pos, m_positions[k], smoothingLength());
    }

    // Boundary contribution
    Eigen::Vector3d boundaryContrib(0.0, 0.0, 0.0);
    for (const BoundarySystem &boundary : boundaries) {
        const size_t boundaryID = boundary.getPointSetID();
        for (size_t j = 0; j < fluidPS.n_neighbors(boundaryID, index); j++) {
            const unsigned int k = fluidPS.neighbor(boundaryID, index, j);
            
            boundaryContrib += m_restDensity * boundary.getParticleVolume(k)
                * particleDensityRatio
                * CubicSpline::gradWeight(pos, boundary.getParticlePos(k), smoothingLength());
        }
    }
    
    return -(fluidContrib + boundaryContrib);
}

Eigen::Vector3d FluidSystem::particleViscosityAcc(size_t index,
                                                  CompactNSearch::NeighborhoodSearch &nsearch,
                                                  const std::vector<BoundarySystem> &boundaries)
{
    const size_t id = getPointSetID();
    const Eigen::Vector3d &pos = getParticlePos(index);
    const double particleMass = this->particleMass();
   
    CompactNSearch::PointSet const& fluidPS = nsearch.point_set(id);
    
    // Fluid contribution
    Eigen::Vector3d fluidContrib (0.0, 0.0, 0.0);
    if (getViscosity() >= 0.0) {
        for (size_t j = 0; j < fluidPS.n_neighbors(id, index); j++) {
            const unsigned int k = fluidPS.neighbor(id, index, j);

            Eigen::Vector3d posDiff = getParticlePos(index) - getParticlePos(k);
            fluidContrib +=
                (particleMass / getParticleDensity(k)) *
                (getParticleVel(index) - getParticleVel(k)) *
                posDiff.dot(CubicSpline::gradWeight(pos, getParticlePos(k), smoothingLength())) /
                (pow(posDiff.norm(), 2) + 0.01 * pow(smoothingLength(), 2));
        }
        fluidContrib *= getViscosity();
    }

    // Boundary contribution
    Eigen::Vector3d boundaryContrib(0.0, 0.0, 0.0);
    for (const BoundarySystem& boundary : boundaries) {
        if (boundary.getViscosity() == 0.0)
            continue;
        
        const size_t boundaryID = boundary.getPointSetID();
        for (size_t j = 0; j < fluidPS.n_neighbors(boundaryID, index); j++) {
            const unsigned int k = fluidPS.neighbor(boundaryID, index, j);

            Eigen::Vector3d posDiff = pos - boundary.getParticlePos(k);
            
            boundaryContrib += boundary.getViscosity() 
                * boundary.getParticleVolume(k) * getParticleVel(index)
                * posDiff.dot(CubicSpline::gradWeight(pos, boundary.getParticlePos(k), smoothingLength()))
                * (pow(posDiff.norm(), 2) + 0.01 * pow(smoothingLength(), 2));
        }
    }
    

    return 2.0 * (fluidContrib + boundaryContrib);
}
