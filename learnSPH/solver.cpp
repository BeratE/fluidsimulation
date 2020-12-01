#include "solver.h"
#include "kernel.h"
#include "util/config.h"
#include "util/vtk_writer.h"
#include <iostream>

#define VECTOR_GRAVITY Eigen::Vector3d(0.0, -9.80665, 0.0)

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
using namespace learnSPH::Kernel;


SolverSPH::SolverSPH(FluidSystem system):
    m_system(system), m_nsearch(1.0)
{    
    m_nsearch.set_radius(Kernel::CubicSpline::support(m_system.smoothingLength()));
    m_system.addToNeighborhood(m_nsearch);
    m_system.initTable();
}

SolverSPH::~SolverSPH()
{
}

double SolverSPH::timeStepCFL()
{
    const double lambda = 0.5;

    auto const &velocities = m_system.getVelocities();
    double maxVelNorm = 0.0;
    for (const auto &vel : velocities) 
        maxVelNorm = maxVelNorm > vel.norm() ? maxVelNorm : vel.norm();
    if (maxVelNorm <= 0.0)
        return m_minTimeStep;
    
    return std::max(m_minTimeStep, lambda * (2*m_system.getParticleRadius() / maxVelNorm));
}

double SolverSPH::integrationStep()
{
    const double deltaT = timeStepCFL();
    
    m_nsearch.find_neighbors();
    m_system.estimateDensity(m_nsearch);
    m_system.updatePressure();
    m_system.updateAcceleration(m_nsearch, m_boundaries);
    
    // Add  gravity
    if (m_gravityEnable) {
        auto &acc = m_system.getAccelerations();
        for (size_t i = 0; i < m_system.getSize(); i++) {
            m_system.setParticleAcc(i, acc[i] + VECTOR_GRAVITY);
        }
    }
    
    semiImplicitEulerStep(deltaT);

    return deltaT;
}

void SolverSPH::addBoundary(BoundarySystem boundary)
{
    m_boundaries.push_back(boundary);
    m_boundaries.back().addToNeighborhood(m_nsearch);
}

void SolverSPH::semiImplicitEulerStep(double deltaT)
{    
    const size_t id = m_system.getPointSetID();
    CompactNSearch::PointSet const& fluidPS = m_nsearch.point_set(id);

    // Update velocities
    const std::vector<Eigen::Vector3d> &velocities = m_system.getVelocities();
    for (size_t i = 0; i < fluidPS.n_points(); i++) {
        Eigen::Vector3d dV = deltaT * m_system.getParticleAcc(i);
        m_system.setParticleVel(i,  velocities[i] + dV);
    }

    // Update positions
    for (size_t i = 0; i < fluidPS.n_points(); i++) {
        const Eigen::Vector3d &fpPos = m_system.getParticlePos(i);
        const Eigen::Vector3d &fpVel = m_system.getParticleVel(i);
        const double fpDensity = m_system.getParticleDensity(i);
        
        Eigen::Vector3d fpVelStar = fpVel;
        if (m_smoothingEnable) { // perform XSPH smoothing
            Eigen::Vector3d fpVelSumOverNeighbors = Eigen::Vector3d::Zero();
            for (size_t j = 0; j < fluidPS.n_neighbors(id, i);
                 j++) {
                const unsigned int k =fluidPS.neighbor(id, i, j);
                fpVelSumOverNeighbors +=
                    (m_system.getParticleVel(k) - fpVel) /
                    (m_system.getParticleDensity(k) + fpDensity) *
                    Kernel::CubicSpline::weight(fpPos, m_system.getParticlePos(k),
                                                m_system.smoothingLength());

                fpVelSumOverNeighbors *= 2 * m_system.particleMass() * m_smoothEps;
                fpVelStar += fpVelSumOverNeighbors;
            }
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
    
    // Single iteration step
    // double deltaT_s = integrationStep();
    // runTime_s += deltaT;
    
    // // Snapshot of first iteration
    // filename.str(std::string());
    // filename <<  SOURCE_DIR << "/res/simulation/" << file << snapShotNr << ".vtk";
    // save_particles_to_vtk(filename.str(), m_system.getPositions(), m_system.getDensities());
    // snapShotNr++;
    // nextSnapShotTime_s = snapShotNr * m_snapShotMS * pow(10, -3);
    // std::cout << "save results to " << filename.str() << std::endl;
    
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
