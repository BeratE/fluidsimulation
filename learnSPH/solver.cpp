#include "solver.h"
#include "kernel.h"
#include "util/config.h"
#include "util/vtk_writer.h"
#include <iostream>

template<typename T>
std::vector<T> interpolateVector(const std::vector<T>& previous,
                                 const std::vector<T>& current,
                                 double prevTime,
                                 double currTime,
                                 double targetTime)
{
    double alpha = (targetTime - prevTime) / (currTime - prevTime);

    auto inter_func =  [alpha](const T& prev, const T& curr)
    {
        return (1.0 - alpha) * prev + alpha * curr;
    };
    
    std::vector<T> interpolation(previous);
    std::transform(previous.begin(), previous.end(),
                   current.begin(), interpolation.begin(), inter_func);
    
    return interpolation;
}


using namespace learnSPH;
using namespace learnSPH::System;
using namespace learnSPH::Kernel;


SolverSPH::SolverSPH(FluidSystem system):
    m_system(system)
{
    mp_nsearch = std::make_shared<CompactNSearch::NeighborhoodSearch>(m_system.getParticleRadius());
    mp_nsearch->set_radius(Kernel::CubicSpline::support(m_system.getSmoothingLength()));
    m_system.addToNeighborhood(mp_nsearch);
}

SolverSPH::~SolverSPH()
{
    mp_nsearch.reset();
}

double SolverSPH::timeStepCFL()
{
    const double lambda = 0.5;

    auto const &velocities = m_system.getVelocities();
    double maxVelNorm = pow(10, -6);
    for (const auto &vel : velocities) 
        maxVelNorm = maxVelNorm > vel.norm() ? maxVelNorm : vel.norm();
    
    return std::min(m_maxTimeStep_s, lambda * (2*m_system.getParticleRadius() / maxVelNorm));
}

double SolverSPH::integrationStep()
{
    const double deltaT = timeStepCFL();
    
    mp_nsearch->find_neighbors();
    
    m_system.updateDensities(m_boundaries);
    m_system.updatePressures(m_param.stiffness);

    applyExternalForces();
    m_system.updateAccelerations(m_boundaries);    
    
    semiImplicitEulerStep(deltaT);

    return deltaT;
}

size_t SolverSPH::addBoundary(const BoundarySystem &boundary)
{
    m_boundaries.push_back(boundary);
    m_boundaries.back().addToNeighborhood(mp_nsearch);
    return m_boundaries.size()-1;
}

void SolverSPH::applyExternalForces()
{
    m_system.clearForces();
    
    Eigen::Vector3d grav_force(0.0, 0.0, 0.0);
    if (m_gravityEnable)
        grav_force = VEC_GRAVITY * m_system.getRestDensity();

    for (size_t i = 0; i < m_system.getSize(); i++) {
        // gravity
        m_system.addParticleForce(i, grav_force);

        // drag force
        auto drag_force = -m_param.drag * m_system.getParticleVel(i);
        m_system.addParticleForce(i, drag_force);        
    }

    // Iterate force objects
    // ...
}

void SolverSPH::semiImplicitEulerStep(double deltaT)
{    
    const size_t id = m_system.getPointSetID();
    CompactNSearch::PointSet const& fluidPS = mp_nsearch->point_set(id);
    
    // Update velocities
    const std::vector<Eigen::Vector3d> &accelerations = m_system.getAccelerations();
    const std::vector<Eigen::Vector3d> &velocities = m_system.getVelocities();
    for (size_t i = 0; i < fluidPS.n_points(); i++) {
        Eigen::Vector3d dV = deltaT * accelerations[i];
        m_system.addParticleVel(i,  dV);
    }

    // Update positions
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
            
            fpVelStar += m_param.smoothing * fpVelSumOverNeighbors;
        }
        m_system.setParticlePos(i, fpPos + deltaT * fpVelStar);
    }
}

void SolverSPH::run(std::string file, double milliseconds)
{
    int boundaryIdx = 0;
    std::stringstream filename;
    for (const BoundarySystem& boundary : m_boundaries) {
        filename.str(std::string());
        filename << SOURCE_DIR << "/res/simulation/"
                 << file << "_boundary" << boundaryIdx++ << ".vtk";
        save_particles_to_vtk(filename.str(), boundary.getPositions(), boundary.getVolumes());
        std::cout << "save results to " << filename.str() << std::endl;
    }
    
    std::vector<Eigen::Vector3d> previousPos(m_system.getPositions());

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

            prevSnapShotTime_s = nextSnapShotTime_s;
            nextSnapShotTime_s = (++snapShotNr) * m_snapShotMS * pow(10, -3);
            
            std::cout << "save results to " << filename.str() << std::endl;
        }
    }
}
