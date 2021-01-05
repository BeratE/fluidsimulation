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


double SolverPBF::integrationStep()
{
    const double deltaT = timeStepCFL();
    
    mp_nsearch->find_neighbors();
    
    m_system.updateDensities(m_boundaries);

    applyExternalForces();
    
    m_system.updateAccelerations(m_boundaries, false, true, true);    
    
    semiImplicitEulerStep(deltaT);
    
    m_system.updateNormalizedDensities();

    return deltaT;
}

void SolverPBF::run(std::string file, double milliseconds, std::vector<Surface::SurfaceInformation>* pOutSurfaceInfos)
{
    std::stringstream filename;
    double runTime_s = 0.0;
    size_t iteration = 0;
    size_t snapShotNr = 0;
    double prevSnapShotTime_s = 0.0;
    double nextSnapShotTime_s = 0.0;
    
    const double END_TIME_s = std::floor(milliseconds / m_snapShotMS)
        * m_snapShotMS * pow(10, -3);
    
    while (runTime_s <= END_TIME_s && ++iteration) {
        std::cout << iteration << " " << runTime_s <<std::endl;
        
        std::vector<Eigen::Vector3d> previousPos(m_system.getPositions());
        
        double deltaT_s = integrationStep();
        runTime_s += deltaT_s;

        mp_nsearch->find_neighbors();

        // constraint solver
        for (size_t k = 0; k < m_npbfIterations; k++) {
            m_system.updateDensities(m_boundaries);

            std::vector<double> lambda;
            for (size_t i = 0; i < m_system.getSize(); i++) {
                lambda.push_back(0.0);
                if (C(i) > 0)
                    lambda[i] = -C(i)/(S(i)+pow(10, -4));
            }
            for (size_t i = 0; i < m_system.getSize(); i++) {                              
                m_system.addParticlePos(i, deltaX(i, lambda));
            }
        }

        // Update velocities
        for (size_t i = 0; i < m_system.getSize(); i++) {
            Eigen::Vector3d v = (m_system.getParticlePos(i) - previousPos[i]) / deltaT_s;
            m_system.setParticleVel(i, v);
        }

        // Snapshot
        if (runTime_s > nextSnapShotTime_s) {
            filename.str(std::string());
            filename <<  SOURCE_DIR << "/res/simulation/" << file << snapShotNr << ".vtk";
            std::vector<Eigen::Vector3d> interpolPos
                = interpolateVector<Eigen::Vector3d>(previousPos,
                                                     m_system.getPositions(),
                                                     prevSnapShotTime_s,
                                                     runTime_s,
                                                     nextSnapShotTime_s);
            save_particles_to_vtk(filename.str(), interpolPos, m_system.getDensities());
           
            std::cout << "save results to " << filename.str() << std::endl;
        }
    }
}

double SolverPBF::C(size_t i)
{
    return (m_system.getParticleDensity(i)/m_system.getRestDensity()) -1;
}

double SolverPBF::S(size_t i)
{
    const size_t id = m_system.getPointSetID();
    const Eigen::Vector3d &pos = m_system.getParticlePos(i);
    const Kernel::CubicSpline::Table& kernelLut = m_system.getKernelLookUp();
    const double inverseVol = m_system.getParticleMass() / m_system.getRestDensity();
    CompactNSearch::PointSet const& fluidPS = mp_nsearch->point_set(id);
    
    double s = 0.0;
    
    // Fluid contribution
    Eigen::Vector3d sum (0.0, 0.0, 0.0);
    sum += inverseVol * kernelLut.gradWeight(pos, pos);
    
    for (size_t j = 0; j < fluidPS.n_neighbors(id, i); j++) {
        const unsigned int k = fluidPS.neighbor(id, i, j);        
        sum += inverseVol * kernelLut.gradWeight(pos, m_system.getParticlePos(k));
    }
    for (const BoundarySystem& boundary : m_boundaries) {
        const size_t boundaryID = boundary.getPointSetID();
        for (size_t j = 0; j < fluidPS.n_neighbors(boundaryID, i); j++) {
            const unsigned int k = fluidPS.neighbor(boundaryID, i, j);
            sum += boundary.getParticleVolume(k) *
                kernelLut.gradWeight(pos, boundary.getParticlePos(k));
        }
    }
    
    s = sum.squaredNorm()/m_system.getParticleMass();
    
    s += (-inverseVol * kernelLut.gradWeight(pos, pos)).squaredNorm()
        / m_system.getParticleMass();
    
    for (size_t j = 0; j < fluidPS.n_neighbors(id, i); j++) {
        const unsigned int k = fluidPS.neighbor(id, i, j);
        s += (-inverseVol * kernelLut.gradWeight(pos, m_system.getParticlePos(k))).squaredNorm()
            / m_system.getParticleMass();
    }
    
    return s;
}

Eigen::Vector3d SolverPBF::deltaX(size_t i, std::vector<double> lambda)
{
    const size_t id = m_system.getPointSetID();
    const Eigen::Vector3d &pos = m_system.getParticlePos(i);
    const Kernel::CubicSpline::Table& kernelLut = m_system.getKernelLookUp();
    CompactNSearch::PointSet const& fluidPS = mp_nsearch->point_set(id);

    Eigen::Vector3d deltaX (0.0, 0.0, 0.0);
    // fluid contribution
    deltaX += 2*lambda[i] * kernelLut.gradWeight(pos, pos);
    
    for (size_t j = 0; j < fluidPS.n_neighbors(id, i); j++) {
        const unsigned int k = fluidPS.neighbor(id, i, j);        
        deltaX += (lambda[i]+lambda[j])
            *kernelLut.gradWeight(pos, m_system.getParticlePos(k));
    }

    deltaX /= m_system.getRestDensity();
    
    // boundary contribution
    for (const BoundarySystem& boundary : m_boundaries) {
        const size_t boundaryID = boundary.getPointSetID();
        for (size_t j = 0; j < fluidPS.n_neighbors(boundaryID, i); j++) {
            const unsigned int k = fluidPS.neighbor(boundaryID, i, j);
            deltaX += (boundary.getParticleVolume(k)/m_system.getParticleMass())
                * lambda[i] * kernelLut.gradWeight(pos, boundary.getParticlePos(k));
        }
    }

    return deltaX;
}
