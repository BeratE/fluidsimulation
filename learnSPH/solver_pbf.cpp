#include "solver_pbf.h"
#include "kernel.h"
#include "config.h"
#include "vtk_writer.h"
#include <iostream>

using namespace learnSPH;
using namespace learnSPH::System;
using namespace learnSPH::Kernel;
using namespace learnSPH::Surface;


SolverPBF::SolverPBF(FluidSystem system):
    Solver(system)
{
}

SolverPBF::~SolverPBF()
{
}

void SolverPBF::integrationStep(double deltaT) {
    mp_nsearch->find_neighbors();
    
    // preparations for the calculations in the semiImplicit Euler integration
    #pragma omp parallel default(none) firstprivate(deltaT)
    {
        m_system.updateDensities(m_boundaries);
        if (m_tensionEnable) {
            m_system.updateNormals();
        }

        updateAccelerations(deltaT);
        semiImplicitEulerStep(deltaT);
    
        #pragma omp single
        {
            // Update estimate based on constraints
            mp_nsearch->find_neighbors();
        }

        for (size_t k = 0; k < m_npbfIterations; k++) {
            #pragma omp barrier
            updatePositionsWithConstraints();           
        }

        // Update Velocities
        #pragma omp for schedule(static)
        for (int i = 0; i < m_system.getSize(); i++) {
            m_system.setParticleVel(
                i, (m_system.getParticlePos(i) - m_system.getParticlePrevPos(i)) / deltaT);
        }
    }
}

void SolverPBF::updateAccelerations(const double deltaT)
{
    initAccelerations();
    
    const size_t id = m_system.getPointSetID();
    CompactNSearch::PointSet const& fluidPS = mp_nsearch->point_set(id);

    // Precalculated density pressure ratios
    static std::vector<double> ratios(m_system.getSize());
    #pragma omp for schedule(static)
    for (int i = 0; i < m_system.getSize(); i++) {
        ratios[i] = (m_system.getParticlePressure(i)
                     / (m_system.getParticleDensity(i)
                        * m_system.getParticleDensity(i)));
    }

    // Synchronization barrier for m_pressureDensityRatios
    #pragma omp barrier 
    
    #pragma omp for schedule(static) 
    for (int i = 0; i < fluidPS.n_points(); i++) {
        // Iterate over fluid neighbors and add contributions to forces 
        for (size_t idx = 0; idx < fluidPS.n_neighbors(id, i); idx++) {
            const unsigned int j = fluidPS.neighbor(id, i, idx);
            updateAccFluidContribution(i, j, ratios[i], ratios[j]);
        }

        // Iterate over boundaries and add contributions to forces
        for (BoundarySystem boundary : m_boundaries) {
            // Iterate over neighboring boundary particles
            const size_t boundaryID = boundary.getPointSetID();
            for (size_t idx = 0; idx < fluidPS.n_neighbors(boundaryID, i); idx++) {
                const unsigned int k = fluidPS.neighbor(boundaryID, i, idx);
                updateAccBoundaryContribution(i, k, ratios[i], boundary);
            }
        }
    }
}


void SolverPBF::updateAccFluidContribution(
    const size_t i,
    const size_t j,
    const double ratio_i,
    const double ratio_j) {

    // add fluid contribution to viscosity acceleration
    m_system.addToParticleAcc(i, m_system.viscAccFluid(i, j));

    // add fluid contribution to tension accelerations. DIVIDE BY PARTICLE MASS
    if (m_tensionEnable)
        m_system.addToParticleAcc(i, m_system.tensionForce(i, j)
                                  / m_system.getParticleMass());
}

void SolverPBF::updateAccBoundaryContribution(
    const size_t i,
    const size_t k,
    const double ratio_i,
    BoundarySystem& boundary) {
    // add boundary contribution to viscosity acceleration
    m_system.addToParticleAcc(i,  m_system.viscAccBoundary(i, k, boundary));
    // add boundary contribution to adhesion accelerations. DIVIDE BY PARTICLE MASS
    if (m_adhesionEnable)
        m_system.addToParticleAcc(i, m_system.adhesionForce(i, k, boundary)
                                  / m_system.getParticleMass());
}


void SolverPBF::updatePositionsWithConstraints()
{
    const size_t fluidID = m_system.getPointSetID();
    CompactNSearch::PointSet const& fluidPS = mp_nsearch->point_set(fluidID);
    
    m_system.updateDensities(m_boundaries);    

    // Compute constraint lambdas
    static std::vector<double> lambdas(m_system.getSize());

    #pragma omp for schedule(static) 
    for (int i = 0; i < m_system.getSize(); i++) {
        lambdas[i] = 0.0;
        if (m_system.getParticleDensity(i) <= m_system.getRestDensity()) {
            continue;
        }                         

        const Eigen::Vector3d pos_i = m_system.getParticlePos(i);
        const double C = m_system.getParticleDensity(i) / m_system.getRestDensity() - 1.0;        
                    
        double S = 0.0;
        Eigen::Vector3d grad = Eigen::Vector3d(0.0, 0.0, 0.0);

        // Fluid neighbors
        for (size_t nIdx = 0; nIdx < fluidPS.n_neighbors(fluidID, i); nIdx++) {
            const unsigned int j = fluidPS.neighbor(fluidID, i, nIdx);
            const Eigen::Vector3d C_grad_j = (m_system.getParticleMass() / m_system.getRestDensity()) 
                * m_system.getKernelLookUp().gradWeight(pos_i, m_system.getParticlePos(j));
            grad += C_grad_j;
            S += C_grad_j.squaredNorm();
        }

        // Boundary neighbors
        for (BoundarySystem boundary : m_boundaries) {
            const size_t boundaryID = boundary.getPointSetID();
            for (size_t bIdx = 0; bIdx < fluidPS.n_neighbors(boundaryID, i); bIdx++) {
                const unsigned int k = fluidPS.neighbor(boundaryID, i, bIdx);
                const Eigen::Vector3d C_grad_k = boundary.getParticleVolume(k)
                    * m_system.getKernelLookUp().gradWeight(pos_i, boundary.getParticlePos(k));
                grad += C_grad_k;
            }
        }
        S += grad.squaredNorm();
        S *= (1.0 / m_system.getParticleMass());
        lambdas[i] = -C / (S + m_eps);
    }

    // Calculate deltaX
    static std::vector<Eigen::Vector3d> deltaX(m_system.getSize());
    
    #pragma omp for schedule(static) 
    for (int i = 0; i < fluidPS.n_points(); i++) {      
        const Eigen::Vector3d pos_i = m_system.getParticlePos(i);
        
        deltaX[i] = 1.0 / m_system.getRestDensity()
            * (lambdas[i] + lambdas[i])
            * m_system.getKernelLookUp().gradWeight(pos_i,
                                                    pos_i);

        
        // Fluid neighbors
        Eigen::Vector3d fluidContrib = Eigen::Vector3d::Zero();
        for (size_t nIdx = 0; nIdx < fluidPS.n_neighbors(fluidID, i); nIdx++) {
            const unsigned int j = fluidPS.neighbor(fluidID, i, nIdx);
            fluidContrib += (lambdas[i] + lambdas[j])
                * m_system.getKernelLookUp().gradWeight(pos_i, m_system.getParticlePos(j));
        }
        deltaX[i] += 1.0 / m_system.getRestDensity() * fluidContrib;
        
        // Boundary neighbors
        for (BoundarySystem boundary : m_boundaries) {
            for (size_t bIdx = 0; bIdx < fluidPS.n_neighbors(boundary.getPointSetID(), i); bIdx++) {
                const unsigned int k = fluidPS.neighbor(boundary.getPointSetID(), i, bIdx);
                deltaX[i] += (boundary.getParticleVolume(k) / m_system.getParticleMass())
                    * lambdas[i] * m_system.getKernelLookUp().gradWeight(pos_i, boundary.getParticlePos(k));
            }
        }
    }

    // Update positions
    #pragma omp for schedule(static) 
    for (int i = 0; i < m_system.getSize(); i++) {
        m_system.addParticlePos(i, deltaX[i]);
    }
}
