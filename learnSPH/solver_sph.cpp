#include "solver_sph.h"
#include "kernel.h"
#include "config.h"
#include "vtk_writer.h"
#include <iostream>
#include <tension/tension.h>

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


double SolverSPH::integrationStep()
{
    const double deltaT = timeStepCFL();
    
    mp_nsearch->find_neighbors();
    
    m_system.updateDensities(m_boundaries);
    m_system.updatePressures(m_stiffness);
    if (m_tensionEnable) {
        m_system.updateNormals();
    }

    applyExternalForces();
    applyTensionForces();
    applyAdhesionForces();
    // drag force
    for (size_t i = 0; i < m_system.getSize(); i++) {
        auto drag_force = -m_drag * m_system.getParticleVel(i);
        m_system.addParticleForce(i, drag_force);        
    }
    
    m_system.updateAccelerations(m_boundaries, true, true, true, true, true);    
    
    semiImplicitEulerStep(deltaT);

    return deltaT;
}

void SolverSPH::run(std::string file, double milliseconds, std::vector<Surface::SurfaceInformation>* pOutSurfaceInfos)
{
    if (pOutSurfaceInfos) {
        mp_nsearch->find_neighbors();
    }
    int boundaryIdx = 0;
    std::stringstream filename;
    for (BoundarySystem& boundary : m_boundaries) {
        filename.str(std::string());
        filename << SOURCE_DIR << "/res/simulation/"
                 << file << "_boundary" << boundaryIdx << ".vtk";
        save_particles_to_vtk(filename.str(), boundary.getPositions(), boundary.getVolumes());
        std::cout << "save results to " << filename.str() << std::endl;
        
        if (pOutSurfaceInfos) {
            boundary.updateNormalizedDensities();
            std::stringstream surfaceFilename;
            surfaceFilename << file << "_boundary_surface" << boundaryIdx;
            pOutSurfaceInfos->push_back(SurfaceInformation(boundary.getPositions(), boundary.getNormalizedDensities(), 
                                        boundary.getKernelLookUp(), boundary.getSmoothingLength(), surfaceFilename.str()));
        }

        boundaryIdx++;
    }


    
    std::vector<Eigen::Vector3d> previousPos(m_system.getPositions());
    std::vector<double> previousNormalizedDensities;
    
    if (pOutSurfaceInfos) {
        m_system.updateNormalizedDensities();
        previousNormalizedDensities = m_system.getNormalizedDensities();
    }


    double runTime_s = 0.0;
    size_t iteration = 0;
    size_t snapShotNr = 0;
    double prevSnapShotTime_s = 0.0;
    double nextSnapShotTime_s = 0.0;
    
    const double END_TIME_s = std::floor(milliseconds / m_snapShotMS)
        * m_snapShotMS * pow(10, -3);
    
    while (runTime_s <= END_TIME_s && ++iteration) {
        std::cout << iteration << " " << runTime_s <<std::endl;
        
        double deltaT_s = integrationStep();
        if (pOutSurfaceInfos) {
            m_system.updateNormalizedDensities();
        }
        runTime_s += deltaT_s;
        
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

            if (pOutSurfaceInfos) {
                std::vector<double> interpolNormalizedDensities
                    = interpolateVector<double>(previousNormalizedDensities,
                    m_system.getNormalizedDensities(),
                    prevSnapShotTime_s,
                    runTime_s,
                    nextSnapShotTime_s);
                std::stringstream surfaceFilename;
                surfaceFilename << file + "_surface" << snapShotNr;
                pOutSurfaceInfos->push_back(SurfaceInformation(
                                                interpolPos, interpolNormalizedDensities,
                                                m_system.getKernelLookUp(),
                                                m_system.getSmoothingLength(),
                                                surfaceFilename.str()));
            }

            prevSnapShotTime_s = nextSnapShotTime_s;
            nextSnapShotTime_s = (++snapShotNr) * m_snapShotMS * pow(10, -3);

            previousPos = m_system.getPositions();
            if (pOutSurfaceInfos) {
                previousNormalizedDensities = m_system.getNormalizedDensities();
            }
            
            std::cout << "save results to " << filename.str() << std::endl;
        }
    }
}
