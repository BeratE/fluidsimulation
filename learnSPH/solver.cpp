#include "solver.h"
#include "kernel.h"
#include "config.h"
#include "vtk_writer.h"
#include <iostream>
#include <surface/surface.h>

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
        grav_force = VEC_GRAVITY * m_system.getRestDensity();

    //#pragma omp parallel for
    for (size_t i = 0; i < m_system.getSize(); i++) {
        // gravity
        m_system.addParticleForce(i, grav_force);     
    }
    // Iterate force objects
    // ...
}

