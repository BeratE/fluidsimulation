#include "solver.h"
#include "kernel.h"
#include "config.h"
#include "vtk_writer.h"
#include <iostream>
#include <surface/surface.h>
#include <tension/tension.h>

using namespace learnSPH;
using namespace learnSPH::System;

Solver::Solver(System::FluidSystem system)
    : m_system(system)
{
    mp_nsearch = std::make_shared<CompactNSearch::NeighborhoodSearch>(
        m_system.getParticleRadius());
    mp_nsearch->set_radius(Kernel::CubicSpline::support(
                               m_system.getSmoothingLength()));
    m_system.addToNeighborhood(mp_nsearch);
}

Solver::~Solver()
{
    mp_nsearch.reset();
}

double Solver::timeStepCFL()
{
    const double lambda = 0.5;

    auto const &velocities = m_system.getVelocities();
    double maxVelNorm = pow(10, -6);
    for (const auto &vel : velocities) 
        maxVelNorm = maxVelNorm > vel.norm() ? maxVelNorm : vel.norm();
    
    return std::min(m_maxTimeStep_s, lambda * (2*m_system.getParticleRadius() / maxVelNorm));
}

size_t Solver::addBoundary(const BoundarySystem &boundary)
{
    m_boundaries.push_back(boundary);
    m_boundaries.back().addToNeighborhood(mp_nsearch);
    return m_boundaries.size()-1;
}

void Solver::applyExternalForces()
{
    m_system.clearForces();
    
    Eigen::Vector3d grav_force(0.0, 0.0, 0.0);
    if (m_gravityEnable)
        grav_force = VEC_GRAVITY * m_system.getParticleMass();

    //#pragma omp parallel for
    for (size_t i = 0; i < m_system.getSize(); i++) {
        // gravity
        m_system.addParticleForce(i, grav_force);     
    }
    // Iterate force objects
    // ...
}

void Solver::applyTensionForces()
{
    using namespace tension;
    if (m_tensionEnable) {
        m_system.clearTensionForces();

        CompactNSearch::PointSet const& fluidPS = mp_nsearch->point_set(m_system.getPointSetID());
        for (size_t i = 0; i < m_system.getSize(); i++) {
            for (size_t j = 0; j < fluidPS.n_neighbors(m_system.getPointSetID(), i); j++) {
                const unsigned int pid = fluidPS.neighbor(m_system.getPointSetID(), i, j);
                m_system.addParticleTensionForce(i, forceTension(
                    m_system.getRestDensity(),
                    m_system.getParticleDensity(i),
                    m_system.getParticleDensity(pid),
                    Cohesion::forceCohesion(
                        m_system.getParticleMass(),
                        m_system.getParticleMass(),
                        m_system.getParticlePos(i),
                        m_system.getParticlePos(pid),
                        m_system.getGamma(),
                        m_system.getC()
                    ),
                    Curvature::forceCurvature(
                        m_system.getParticleMass(),
                        m_system.getParticleNormal(i),
                        m_system.getParticleNormal(pid),
                        m_system.getC()
                    )
                ));
            }
        }
    }
}

void Solver::applyAdhesionForces()
{
    using namespace tension;
    if (m_adhesionEnable) {
        m_system.clearAdhesionForces();

        CompactNSearch::PointSet const& fluidPS = mp_nsearch->point_set(m_system.getPointSetID());
        for (size_t i = 0; i < m_system.getSize(); i++) {
            for (System::BoundarySystem boundary : m_boundaries) {
                for (size_t k = 0; k < fluidPS.n_neighbors(boundary.getPointSetID(), i); k++) {
                    const unsigned int pid = fluidPS.neighbor(boundary.getPointSetID(), i, k);
                    m_system.addParticleAdhesionForce(i, Adhesion::forceAdhesion(
                        m_system.getParticleMass(),
                        boundary.getParticleVolume(pid),
                        m_system.getParticlePos(i),
                        boundary.getParticlePos(pid),
                        boundary.getBeta(),
                        m_system.getC()
                    ));
                }
            }
        }
    }
}

void Solver::semiImplicitEulerStep(double deltaT)
{    
    const size_t id = m_system.getPointSetID();
    CompactNSearch::PointSet const& fluidPS = mp_nsearch->point_set(id);
    
    // Update velocities
    const std::vector<Eigen::Vector3d> &accelerations = m_system.getAccelerations();
    const std::vector<Eigen::Vector3d> &velocities = m_system.getVelocities();
    
    //#pragma omp parallel for
    for (size_t i = 0; i < fluidPS.n_points(); i++) {
        Eigen::Vector3d dV = deltaT * accelerations[i];
        m_system.addParticleVel(i,  dV);
    }

    // Update positions
    //#pragma omp parallel for
    for (size_t i = 0; i < fluidPS.n_points(); i++) {
        const Eigen::Vector3d &fpPos = m_system.getParticlePos(i);
        const Eigen::Vector3d &fpVel = m_system.getParticleVel(i);
        const double fpDensity = m_system.getParticleDensity(i);        
        Eigen::Vector3d fpVelStar = fpVel;
        
        // perform XSPH smoothing
        if (m_smoothingEnable) { 
            Eigen::Vector3d fpVelSumOverNeighbors(0.0, 0.0, 0.0);
            
            for (size_t j = 0; j < fluidPS.n_neighbors(id, i); j++) {
                const unsigned int k = fluidPS.neighbor(id, i, j);
                fpVelSumOverNeighbors +=
                    (2 * m_system.getParticleMass() *
                    (m_system.getParticleVel(k) - fpVel)  *
                    Kernel::CubicSpline::weight(fpPos, m_system.getParticlePos(k),
                                                m_system.getSmoothingLength()))
                    / (m_system.getParticleDensity(k) + fpDensity);
            }
            
            fpVelStar += m_xsphSmoothing * fpVelSumOverNeighbors;
        }
        m_system.setParticlePos(i, fpPos + deltaT * fpVelStar);
    }
}

