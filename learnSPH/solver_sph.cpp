#include "solver_sph.h"
#include "kernel.h"
#include "config.h"
#include "vtk_writer.h"
#include <iostream>

using namespace learnSPH;
using namespace learnSPH::System;
using namespace learnSPH::Kernel;
using namespace learnSPH::Surface;

SolverSPH::SolverSPH(FluidSystem system):
    Solver(system)
{
}

SolverSPH::~SolverSPH()
{
}

double SolverSPH::integrationStep(const std::vector<Eigen::Vector3d>& previousPos)
{
    const double deltaT = timeStepCFL();
    
    mp_nsearch->find_neighbors();
    
    // preparations for the calculations in the semiImplicit Euler integration
    m_system.updateDensities(m_boundaries);
    m_system.updatePressures(m_stiffness);
    if (m_tensionEnable) {
        m_system.updateNormals();
    }

    updateAccelerations(deltaT);
    semiImplicitEulerStep(deltaT);

    return deltaT;
}


void SolverSPH::updateAccelerations(double deltaT)
{
    initAccelerations();
    
    //densities, pressures and normals are already calculated 
    const size_t id = m_system.getPointSetID();
    CompactNSearch::PointSet const& fluidPS = mp_nsearch->point_set(id);

    // Used in calculation of pressureAccelerations
    std::vector<double> pressureDensityRatios;
    for (size_t i = 0; i < m_system.getSize(); i++) {
        pressureDensityRatios.push_back(m_system.getParticlePressure(i)
                                        / (m_system.getParticleDensity(i)
                                           * m_system.getParticleDensity(i)));
    }
    
    #pragma omp parallel for schedule(static) 
    for (int i = 0; i < fluidPS.n_points(); i++) {
        // Contributions of pair (i, i) - NOT NECESSARY FOR VISCOSITY AND TENSION
        m_system.addToParticleAcc(i, -m_system.pressureAccFluid(i, i,
                                                                pressureDensityRatios[i],
                                                                pressureDensityRatios[i]));

        // Iterate over fluid neighbors and add contributions to forces 
        for (size_t idx = 0; idx < fluidPS.n_neighbors(id, i); idx++) {
            const unsigned int j = fluidPS.neighbor(id, i, idx);
            updateAccFluidContribution(i, j, pressureDensityRatios[i],
                                       pressureDensityRatios[j]);
        }

        // Iterate over boundaries and add contributions to forces
        for (BoundarySystem boundary : m_boundaries) {
            // Iterate over neighboring boundary particles
            const size_t boundaryID = boundary.getPointSetID();
            for (size_t idx = 0; idx < fluidPS.n_neighbors(boundaryID, i); idx++) {
                const unsigned int k = fluidPS.neighbor(boundaryID, i, idx);
                updateAccBoundaryContribution(i, k, pressureDensityRatios[i], boundary);
            }
        }
    }
}


void SolverSPH::updateAccFluidContribution(const size_t i, 
                                           const size_t j, 
                                           const double ratio_i, 
                                           const double ratio_j) {
    // add fluid contribution to pressure acceleration
    m_system.addToParticleAcc(i, -m_system.pressureAccFluid(i, j,
                                                         ratio_i,
                                                         ratio_j));

    // add fluid contribution to viscosity acceleration
    m_system.addToParticleAcc(i, m_system.viscAccFluid(i, j));

    // add fluid contribution to tension accelerations. DIVIDE BY PARTICLE MASS
    if (m_tensionEnable)
        m_system.addToParticleAcc(i, m_system.tensionForce(i, j)
                                  / m_system.getParticleMass());
}

void SolverSPH::updateAccBoundaryContribution(const size_t i,
                                              const size_t k,
                                              const double ratio_i,
                                              BoundarySystem& boundary) {
    // add boundary contribution to pressure acceleration
    m_system.addToParticleAcc(i, -m_system.pressureAccBoundary(i, k, ratio_i, boundary));
    // add boundary contribution to viscosity acceleration
     m_system.addToParticleAcc(i, m_system.viscAccBoundary(i, k, boundary));
    // add boundary contribution to adhesion accelerations. DIVIDE BY PARTICLE MASS
    if (m_adhesionEnable)
        m_system.addToParticleAcc(i, m_system.adhesionForce(i, k, boundary) / m_system.getParticleMass());
}
