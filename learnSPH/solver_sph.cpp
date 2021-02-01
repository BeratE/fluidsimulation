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


void SolverSPH::newRun(std::string file, double milliseconds, std::vector<Surface::SurfaceInformation>* pOutSurfaceInfos)
{
    // Store information regarding boundaries
    int boundaryIdx = 0;
    std::stringstream filename;
    for (BoundarySystem& boundary : m_boundaries) {
        // Write boundary particles to file
        filename.str(std::string());
        filename << SOURCE_DIR << "/res/simulation/"
            << file << "_boundary" << boundaryIdx << ".vtk";
        save_particles_to_vtk(filename.str(), boundary.getPositions(), boundary.getVolumes());
        std::cout << "save results to " << filename.str() << std::endl;

        boundaryIdx++;
    }

    // Required for interpolation between simulation steps 
    std::vector<Eigen::Vector3d> previousPos(m_system.getPositions());

    // Controll variables
    double runTime_s = 0.0;
    size_t iteration = 0;
    size_t snapShotNr = 0;
    double prevTime_s = 0.0;
    double nextSnapShotTime_s = 0.0;
    const double END_TIME_s = std::floor(milliseconds / m_snapShotMS) * m_snapShotMS * pow(10, -3);

    while (runTime_s <= END_TIME_s && ++iteration) {
        std::cout << iteration << " " << runTime_s << std::endl;
        double deltaT_s = newIntegrationStep();
        runTime_s += deltaT_s;

        if (runTime_s > nextSnapShotTime_s) {
            filename.str(std::string());
            filename << SOURCE_DIR << "/res/simulation/" << file << snapShotNr << ".vtk";
            std::vector<Eigen::Vector3d> interpolPos
                = interpolateVector<Eigen::Vector3d>(previousPos,
                    m_system.getPositions(),
                    prevTime_s,
                    runTime_s,
                    nextSnapShotTime_s);
            save_particles_to_vtk(filename.str(), interpolPos, m_system.getDensities());

            nextSnapShotTime_s = (++snapShotNr) * m_snapShotMS * pow(10, -3);

            std::cout << "save results to " << filename.str() << std::endl;
        }
            prevTime_s = runTime_s;
            previousPos = m_system.getPositions();
    }
}

void SolverSPH::newSemiImplicitEulerStep(double deltaT)
{
    //densities, pressures and normals are already calculated 
    const Eigen::Vector3d GRAVITATIONAL_ACC(0, -9.80665, 0);

    const size_t id = m_system.getPointSetID();
    CompactNSearch::PointSet const& fluidPS = mp_nsearch->point_set(id);

    // Intermediate results
    std::vector<Eigen::Vector3d> accelerations;

    // Used in calculation of pressureAccelerations
    std::vector<double> pressureDensityRatios;
    for (size_t i = 0; i < m_system.getSize(); i++) {
        pressureDensityRatios.push_back(m_system.getParticlePressure(i) / (m_system.getParticleDensity(i) * m_system.getParticleDensity(i)));
    }

    ////////////////////////////////////////////////////////////////////////
    // Calculate Accelerations
    ////////////////////////////////////////////////////////////////////////

    // If gravity is enabled, we can add gravitational acceleration to every particle
    for (int i = 0; i < fluidPS.n_points(); i++) {
        if (m_gravityEnable) {
            accelerations.push_back(GRAVITATIONAL_ACC);
        }
        else {
            accelerations.push_back(Eigen::Vector3d::Zero());
        }
    }
    #pragma omp parallel for schedule(static) num_threads(omp_get_num_procs())
    for (int i = 0; i < fluidPS.n_points(); i++) {
        // Contributions of pair (i, i) - NOT NECESSARY FOR VISCOSITY AND TENSION
        accelerations[i] -= m_system.pressureAccFluid(i, i,
            pressureDensityRatios[i],
            pressureDensityRatios[i]);

        // Iterate over fluid neighbors and add contributions to forces 
        for (size_t idx = 0; idx < fluidPS.n_neighbors(m_system.getPointSetID(), i); idx++) {
            const unsigned int j = fluidPS.neighbor(m_system.getPointSetID(), i, idx);
            updateAccFluidContribution(accelerations, i, j, pressureDensityRatios[i], pressureDensityRatios[j]);
        }

        // Iterate over boundaries and add contributions to forces
        for (BoundarySystem boundary : m_boundaries) {
            // Iterate over neighboring boundary particles
            for (size_t idx = 0; idx < fluidPS.n_neighbors(boundary.getPointSetID(), i); idx++) {
                const unsigned int k = fluidPS.neighbor(boundary.getPointSetID(), i, idx);
                updateAccBoundaryContribution(accelerations, i, k, pressureDensityRatios[i], boundary);
            }
        }
    }

    ////////////////////////////////////////////////////////////////////////
    // Update Velocities
    ////////////////////////////////////////////////////////////////////////
    #pragma omp parallel for schedule(static) num_threads(omp_get_num_procs())
    for (int i = 0; i < m_system.getSize(); i++) {
        m_system.updateVelocity(i, m_system.getParticleVel(i) + deltaT * accelerations[i]);
    }

    ////////////////////////////////////////////////////////////////////////
    // Update Positions
    ////////////////////////////////////////////////////////////////////////
    #pragma omp parallel for schedule(static) num_threads(omp_get_num_procs())
    for (int i = 0; i < fluidPS.n_points(); i++) {
        if (m_smoothingEnable) {
            // Calculate smoothingTerm from neighborhood
            Eigen::Vector3d smoothingTerm = Eigen::Vector3d::Zero();
            for (size_t idx = 0; idx < fluidPS.n_neighbors(m_system.getPointSetID(), i); idx++) {
                const unsigned int j = fluidPS.neighbor(m_system.getPointSetID(), i, idx);
                smoothingTerm += m_system.smoothingTerm(i, j);
            }
            m_system.updatePosition(i, m_system.getParticlePos(i) + deltaT * (m_system.getParticleVel(i) + m_xsphSmoothing * smoothingTerm));
        }
        else {
            m_system.updatePosition(i, m_system.getParticlePos(i) + deltaT * m_system.getParticleVel(i));
        }
    }

}

double SolverSPH::newIntegrationStep()
{
    const double deltaT = timeStepCFL();
    
    mp_nsearch->find_neighbors();
    
    // preparations for the calculations in the semiImplicit Euler integration
    m_system.updateDensities(m_boundaries);
    m_system.updatePressures(m_stiffness);
    if (m_tensionEnable) {
        m_system.updateNormals();
    }
    // TODO: Maybe add drag forces
    newSemiImplicitEulerStep(deltaT);

    return deltaT;
}

void SolverSPH::updateAccFluidContribution(std::vector<Eigen::Vector3d>& accelerations, 
    const size_t i, 
    const size_t j, 
    const double ratio_i, 
    const double ratio_j) {
    // add fluid contribution to pressure acceleration
    accelerations[i] -= m_system.pressureAccFluid(i, j,
        ratio_i,
        ratio_j);

    // add fluid contribution to viscosity acceleration
    accelerations[i] += m_system.viscAccFluid(i, j);

    // add fluid contribution to tension accelerations. DIVIDE BY PARTICLE MASS
    if (m_tensionEnable)
        accelerations[i] += m_system.tensionForce(i, j) / m_system.getParticleMass();
}

void SolverSPH::updateAccBoundaryContribution(std::vector<Eigen::Vector3d>& accelerations,
    const size_t i,
    const size_t k,
    const double ratio_i,
    BoundarySystem& boundary) {
    // add boundary contribution to pressure acceleration
    accelerations[i] -= m_system.pressureAccBoundary(i, k, ratio_i, boundary);
    // add boundary contribution to viscosity acceleration
    accelerations[i] += m_system.viscAccBoundary(i, k, boundary);
    // add boundary contribution to adhesion accelerations. DIVIDE BY PARTICLE MASS
    if (m_adhesionEnable)
        accelerations[i] += m_system.adhesionForce(i, k, boundary) / m_system.getParticleMass();
}




