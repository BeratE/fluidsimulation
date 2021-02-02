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

void SolverPBF::newRun(std::string file, double milliseconds,
    std::vector<Surface::SurfaceInformation>* pOutSurfaceInfos) {

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
        double deltaT_s = newIntegrationStep(previousPos);
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

double SolverPBF::newIntegrationStep(std::vector<Eigen::Vector3d>& previousPos) {
    const double deltaT = timeStepCFL();

    mp_nsearch->find_neighbors();

    // preparations for the calculations in the semiImplicit Euler integration
    m_system.updateDensities(m_boundaries);
    if (m_tensionEnable) {
        m_system.updateNormals();
    }

    newSemiImplicitEulerStep(deltaT);

    // Update estimate based on constraints
    mp_nsearch->find_neighbors(); // TODO: Maybe do for every iteration
    
    for (size_t iterationNr = 0; iterationNr < m_npbfIterations; iterationNr++) {
        updatePositionsWithConstraints();
    }

    // Update Velocities
    #pragma omp parallel for schedule(static) num_threads(omp_get_num_procs())
    for (int i = 0; i < m_system.getSize(); i++) {
        m_system.updateVelocity(i,(m_system.getParticlePos(i) - previousPos[i]) / deltaT);
    }
        

    return deltaT;
}

void SolverPBF::newSemiImplicitEulerStep(const double deltaT) {
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
    for (size_t i = 0; i < fluidPS.n_points(); i++) {
        if (m_gravityEnable) {
            accelerations.push_back(GRAVITATIONAL_ACC);
        }
        else {
            accelerations.push_back(Eigen::Vector3d::Zero());
        }
    }
    #pragma omp parallel for schedule(static) num_threads(omp_get_num_procs())
    for (int i = 0; i < fluidPS.n_points(); i++) {
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
    
    if (m_smoothingEnable) {
        std::vector<Eigen::Vector3d> smoothingTerms(m_system.getSize(), Eigen::Vector3d::Zero());
        #pragma omp parallel for schedule(static) num_threads(omp_get_num_procs())
        for (int i = 0; i < fluidPS.n_points(); i++) {
            for (size_t idx = 0; idx < fluidPS.n_neighbors(m_system.getPointSetID(), i); idx++) {
                const unsigned int j = fluidPS.neighbor(m_system.getPointSetID(), i, idx);
                smoothingTerms[i] += m_system.smoothingTerm(i, j);
            }
        }
        #pragma omp parallel for schedule(static) num_threads(omp_get_num_procs())
        for (int i = 0; i < fluidPS.n_points(); i++) {
            m_system.updatePosition(i, m_system.getParticlePos(i) + deltaT * (m_system.getParticleVel(i) + m_xsphSmoothing * smoothingTerms[i]));
        }
    }
    else {
        #pragma omp parallel for schedule(static) num_threads(omp_get_num_procs())
        for (int i = 0; i < fluidPS.n_points(); i++) {
            m_system.updatePosition(i, m_system.getParticlePos(i) + deltaT * m_system.getParticleVel(i));
        }
    }
}

void SolverPBF::updateAccFluidContribution(std::vector<Eigen::Vector3d>& accelerations,
    const size_t i,
    const size_t j,
    const double ratio_i,
    const double ratio_j) {

    // add fluid contribution to viscosity acceleration
    accelerations[i] += m_system.viscAccFluid(i, j);

    // add fluid contribution to tension accelerations. DIVIDE BY PARTICLE MASS
    if (m_tensionEnable)
        accelerations[i] += m_system.tensionForce(i, j) / m_system.getParticleMass();
}

void SolverPBF::updateAccBoundaryContribution(std::vector<Eigen::Vector3d>& accelerations,
    const size_t i,
    const size_t k,
    const double ratio_i,
    BoundarySystem& boundary) {
    // add boundary contribution to viscosity acceleration
    accelerations[i] += m_system.viscAccBoundary(i, k, boundary);
    // add boundary contribution to adhesion accelerations. DIVIDE BY PARTICLE MASS
    if (m_adhesionEnable)
        accelerations[i] += m_system.adhesionForce(i, k, boundary) / m_system.getParticleMass();
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
    {
        #pragma omp parallel for schedule(static) num_threads(omp_get_num_procs())
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
                    const Eigen::Vector3d C_grad_k = boundary.getParticleVolume(k) * m_system.getKernelLookUp().gradWeight(pos_i, boundary.getParticlePos(k));
                    grad += C_grad_k;
                }
            }
            S += 1.0 / m_system.getParticleMass() * grad.norm() * grad.norm();
            lambdas[i] = -C / (S + m_eps);
        }
    }
    // Calculate deltaX
    std::vector<Eigen::Vector3d> deltaX(m_system.getSize(), Eigen::Vector3d::Zero());
    CompactNSearch::PointSet const& fluidPS = mp_nsearch->point_set(m_system.getPointSetID());
    #pragma omp parallel for schedule(static) num_threads(omp_get_num_procs())
    for (int i = 0; i < fluidPS.n_points(); i++) {
        // Fluid neighbors
        Eigen::Vector3d fluidContrib = Eigen::Vector3d::Zero();
        Eigen::Vector3d pos_i = m_system.getParticlePos(i);
        for (size_t nIdx = 0; nIdx < fluidPS.n_neighbors(m_system.getPointSetID(), i); nIdx++) {
            const unsigned int j = fluidPS.neighbor(m_system.getPointSetID(), i, nIdx);
            fluidContrib += ( ( m_system.getParticleMass() / m_system.getParticleMass() ) * lambdas[i] + lambdas[j]) * m_system.getKernelLookUp().gradWeight(pos_i, m_system.getParticlePos(j));
        }
        deltaX[i] += 1.0 / m_system.getRestDensity() * fluidContrib;
        // Boundary neighbors
        for (BoundarySystem boundary : m_boundaries) {
            for (size_t bIdx = 0; bIdx < fluidPS.n_neighbors(boundary.getPointSetID(), i); bIdx++) {
                const unsigned int k = fluidPS.neighbor(boundary.getPointSetID(), i, bIdx);
                deltaX[i] += (boundary.getParticleVolume(k) / m_system.getParticleMass()) * lambdas[i] * m_system.getKernelLookUp().gradWeight(pos_i, boundary.getParticlePos(k));
            }
        }
    }

    // Update positions
    #pragma omp parallel for schedule(static) num_threads(omp_get_num_procs())
    for (int i = 0; i < m_system.getSize(); i++) {
        m_system.addParticlePos(i, deltaX[i]);
    }
}
