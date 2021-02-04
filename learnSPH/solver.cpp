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
    double maxVelNorm = 0.000001; // some small number
    for (const auto &vel : velocities) 
        if (maxVelNorm < vel.norm()) {
            maxVelNorm = vel.norm();
        }
    return std::min(m_maxTimeStep_s, lambda * (2*m_system.getParticleRadius() / maxVelNorm));
}

size_t Solver::addBoundary(const BoundarySystem &boundary)
{
    m_boundaries.push_back(boundary);
    m_boundaries.back().addToNeighborhood(mp_nsearch);
    return m_boundaries.size()-1;
}


void Solver::run(
    std::string file, double milliseconds,
    std::vector<Surface::SurfaceInformation>* pOutSurfaceInfos) {

    if (pOutSurfaceInfos) {
        mp_nsearch->find_neighbors();
    }
    
    // Store information regarding boundaries
    int boundaryIdx = 0;
    std::stringstream filename;
    std::stringstream surfaceFilename;
        
    for (BoundarySystem& boundary : m_boundaries) {
        // Write boundary particles to file
        filename.str(std::string());
        filename << file << "_boundary" << boundaryIdx << ".vtk";
        save_particles_to_vtk(filename.str(),
                              boundary.getPositions(), boundary.getVolumes());

         if (pOutSurfaceInfos) {
             boundary.updateNormalizedDensities();
             surfaceFilename.str(std::string());
             surfaceFilename << file << "_boundary_surface" << boundaryIdx;
             pOutSurfaceInfos->push_back(
                 Surface::SurfaceInformation(
                     boundary.getPositions(), boundary.getNormalizedDensities(), 
                     boundary.getKernelLookUp(), boundary.getSmoothingLength(),
                     surfaceFilename.str()));
         }

        boundaryIdx++;
    }
    

    // Required for interpolation between simulation steps 
    std::vector<Eigen::Vector3d> previousPos(m_system.getPositions());
    std::vector<double> previousNormalizedDensities;

    if (pOutSurfaceInfos) {
        m_system.updateNormalizedDensities();
        previousNormalizedDensities = m_system.getNormalizedDensities();
    }

    // Controll variables
    double runTime_s = 0.0;
    size_t iteration = 0;
    size_t snapShotNr = 0;
    double prevTime_s = 0.0;
    double nextSnapShotTime_s = 0.0;
    const double END_TIME_s = std::floor(milliseconds / m_snapShotMS)
        * m_snapShotMS * pow(10, -3);

    
    while (runTime_s <= END_TIME_s && ++iteration) {
        std::cout << iteration << " " << runTime_s << std::endl;
        // Propagate System
        double deltaT_s = timeStepCFL();
        integrationStep(deltaT_s, previousPos);
        runTime_s += deltaT_s;

        // Take Snapshot
        if (runTime_s > nextSnapShotTime_s) {
            filename.str(std::string());
            filename << file << snapShotNr << ".vtk";
            std::vector<Eigen::Vector3d> interpolPos
                = interpolateVector<Eigen::Vector3d>(previousPos,
                                                     m_system.getPositions(),
                                                     prevTime_s,
                                                     runTime_s,
                                                     nextSnapShotTime_s);
            save_particles_to_vtk(filename.str(), interpolPos, m_system.getDensities());

            if (pOutSurfaceInfos) {
                surfaceFilename.str(std::string());
                surfaceFilename << file + "_surface" << snapShotNr;
                std::vector<double> interpolNormalizedDensities
                    = interpolateVector<double>(previousNormalizedDensities,
                                                m_system.getNormalizedDensities(),
                                                prevTime_s,
                                                runTime_s,
                                                nextSnapShotTime_s);                
                pOutSurfaceInfos->push_back(Surface::SurfaceInformation(
                                                interpolPos, interpolNormalizedDensities,
                                                m_system.getKernelLookUp(),
                                                m_system.getSmoothingLength(),
                                                surfaceFilename.str()));
            }

            nextSnapShotTime_s = (++snapShotNr) * m_snapShotMS*0.001;            
        }
        prevTime_s = runTime_s;
        previousPos = m_system.getPositions();
        if (pOutSurfaceInfos) {
            previousNormalizedDensities = m_system.getNormalizedDensities();
        }
    }
}

void Solver::semiImplicitEulerStep(double deltaT) {
    const size_t id = m_system.getPointSetID();
    CompactNSearch::PointSet const& fluidPS = mp_nsearch->point_set(id);
    
    ////////////////////////////////////////////////////////////////////////
    // Update Velocities
    ////////////////////////////////////////////////////////////////////////
    #pragma omp for schedule(static) 
    for (int i = 0; i < m_system.getSize(); i++) {
        m_system.addToParticleVel(i, deltaT * m_system.getParticleAcc(i));
    }

    ////////////////////////////////////////////////////////////////////////
    // Update Positions
    ////////////////////////////////////////////////////////////////////////    
    if (m_smoothingEnable) {
        static std::vector<Eigen::Vector3d> smoothingTerms(fluidPS.n_points());

        #pragma omp for schedule(static)
        for (int i = 0; i < fluidPS.n_points(); i++) {
            smoothingTerms[i] = Eigen::Vector3d(0.0, 0.0, 0.0);
            for (size_t idx = 0; idx < fluidPS.n_neighbors(m_system.getPointSetID(), i); idx++) {
                const unsigned int j = fluidPS.neighbor(m_system.getPointSetID(), i, idx);
                smoothingTerms[i] += m_system.smoothingTerm(i, j);
            }
        }
 
        #pragma omp barrier // synchronization of smoothingTerms
            
        #pragma omp for schedule(static)
        for (int i = 0; i < fluidPS.n_points(); i++) {
            m_system.addToParticlePos(
                i, deltaT * (m_system.getParticleVel(i) +
                             m_xsphSmoothing * smoothingTerms[i]));
        }

    }
    else {
        #pragma omp for schedule(static)
        for (int i = 0; i < fluidPS.n_points(); i++) {
            m_system.addToParticlePos(i, deltaT * m_system.getParticleVel(i));
        }
    }
}

void Solver::initAccelerations()
{
    // If gravity is enabled, we can add gravitational acceleration to every particle
    #pragma omp for schedule(static) 
    for (int i = 0; i < m_system.getSize(); i++) {
        Eigen::Vector3d acc(Eigen::Vector3d::Zero());
        if (m_gravityEnable) {
            acc += VEC_GRAVITY;
        }

        m_system.setParticleAcc(i, acc);
    }
}
