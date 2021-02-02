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

double SolverPBF::integrationStep(const std::vector<Eigen::Vector3d>& previousPos) {
    const double deltaT = timeStepCFL();

    mp_nsearch->find_neighbors();

    // preparations for the calculations in the semiImplicit Euler integration
    m_system.updateDensities(m_boundaries);
    if (m_tensionEnable) {
        m_system.updateNormals();
    }

    updateAccelerations(deltaT);
    semiImplicitEulerStep(deltaT);

    // Update estimate based on constraints
    mp_nsearch->find_neighbors(); // TODO: Maybe do for every iteration
    
    for (size_t iterationNr = 0; iterationNr < m_npbfIterations; iterationNr++) {
        updatePositionsWithConstraints();
    }

    // Update Velocities
    #pragma omp parallel for schedule(static) 
    for (int i = 0; i < m_system.getSize(); i++) {
        m_system.setParticleVel(i, (m_system.getParticlePos(i) - previousPos[i]) / deltaT);
    }        

    return deltaT;
}

void SolverPBF::updateAccelerations(const double deltaT) {
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

void SolverPBF::updatePositionsWithConstraints() {
    m_system.updateDensities(m_boundaries);
    std::vector<int> high_density_particles;
    high_density_particles.reserve(m_system.getSize());
    for (size_t i = 0; i < m_system.getSize(); i++) {
        if (m_system.getParticleDensity(i) > m_system.getRestDensity()) {
            high_density_particles.push_back(i);
        }
    }

    // Compute constraint lambdas
    std::vector<double> lambdas(m_system.getSize(), 0.0);

    #pragma omp parallel for schedule(static) 
    for (int idx = 0; idx < high_density_particles.size(); idx++) {
        const size_t i = high_density_particles[idx];
        const Eigen::Vector3d pos_i = m_system.getParticlePos(i);
        const double C = m_system.getParticleDensity(i) / m_system.getRestDensity() - 1.0;
        Eigen::Vector3d grad = Eigen::Vector3d::Zero();
        double S = 0.0;

            
        CompactNSearch::PointSet const& fluidPS = mp_nsearch->point_set(m_system.getPointSetID());
            
        // Fluid neighbors
        for (size_t nIdx = 0; nIdx < fluidPS.n_neighbors(m_system.getPointSetID(), i); nIdx++) {
            const unsigned int j = fluidPS.neighbor(m_system.getPointSetID(), i, nIdx);
            const Eigen::Vector3d C_grad_j = ( m_system.getParticleMass() / m_system.getRestDensity() ) 
                * m_system.getKernelLookUp().gradWeight(pos_i, m_system.getParticlePos(j));
            grad += C_grad_j;
            S += 1.0 / m_system.getParticleMass() * C_grad_j.norm() * C_grad_j.norm();
        }

        // Boundary neighbors
        for (BoundarySystem boundary : m_boundaries) {
            for (size_t bIdx = 0; bIdx < fluidPS.n_neighbors(boundary.getPointSetID(), i); bIdx++) {
                const unsigned int k = fluidPS.neighbor(boundary.getPointSetID(), i, bIdx);
                const Eigen::Vector3d C_grad_k = boundary.getParticleVolume(k)
                    * m_system.getKernelLookUp().gradWeight(pos_i, boundary.getParticlePos(k));
                grad += C_grad_k;
            }
        }
        S += 1.0 / m_system.getParticleMass() * grad.norm() * grad.norm();
        lambdas[i] = -C / (S + m_eps);
    }

    // Calculate deltaX
    std::vector<Eigen::Vector3d> deltaX(m_system.getSize(), Eigen::Vector3d::Zero());
    CompactNSearch::PointSet const& fluidPS = mp_nsearch->point_set(m_system.getPointSetID());
    
    #pragma omp parallel for schedule(static) 
    for (int i = 0; i < fluidPS.n_points(); i++) {
        // Fluid neighbors
        Eigen::Vector3d fluidContrib = Eigen::Vector3d::Zero();
        Eigen::Vector3d pos_i = m_system.getParticlePos(i);
        for (size_t nIdx = 0; nIdx < fluidPS.n_neighbors(m_system.getPointSetID(), i); nIdx++) {
            const unsigned int j = fluidPS.neighbor(m_system.getPointSetID(), i, nIdx);
            fluidContrib += (( m_system.getParticleMass() / m_system.getParticleMass())
                             * lambdas[i] + lambdas[j])
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
    #pragma omp parallel for schedule(static) 
    for (int i = 0; i < m_system.getSize(); i++) {
        m_system.addParticlePos(i, deltaX[i]);
    }
}
