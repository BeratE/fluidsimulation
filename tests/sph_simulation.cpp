#include "catch.hpp"
#include "config.h"
#include "vtk_writer.h"
#include "learnSPH/kernel.h"
#include "learnSPH/solver_sph.h"
#include "learnSPH/system/emitter.h"
#include "learnSPH/system/particlesystem.h"

using namespace learnSPH;
using namespace learnSPH::System;

TEST_CASE("semiEuler", "Test semi-implicit Euler for a cube of particles with gravity only") {
    std::cout << "Testing semi-implicit Euler.." << std::endl;

    const size_t N_STEPS = 100;
    const double particleDiameter = 0.1;

    // Sample Particles in a Box
    FluidSystem particles = Emitter().sampleFluidBox(Eigen::Vector3d(0, 0, 0),
                                                     Eigen::Vector3d(1, 1, 1),
                                                     particleDiameter);

    std::stringstream filename;
    SolverSPH solver(particles);
    solver.enableGravity(false);
    
    SECTION("EulerSmoothON") {
        solver.enableSmoothing(true);
        
        filename << SOURCE_DIR << "/res/integration/test_semi_euler_smooth_on"
                 << 0 << ".vtk";
        save_particles_to_vtk(filename.str(), particles.getPositions(),
                              solver.getSystem().getDensities());

        for (size_t steps = 1; steps < N_STEPS; steps++) {
            filename.str("");
            filename << SOURCE_DIR << "/res/integration/test_semi_euler_smooth_on"
                     << steps << ".vtk";

            solver.integrationStep();
            
            save_particles_to_vtk(filename.str(), solver.getSystem().getPositions(),
                                  solver.getSystem().getDensities());

            std::cout << "Results saved to " << filename.str() << std::endl;
        }
    }
    
    SECTION("EulerSmoothOFF") {
        
        solver.enableSmoothing(false);
        
        filename << SOURCE_DIR << "/res/integration/test_semi_euler_smooth_on"
                 << 0 << ".vtk";
        save_particles_to_vtk(filename.str(), particles.getPositions(),
                              solver.getSystem().getDensities());

        for (size_t steps = 1; steps < N_STEPS; steps++) {
            filename.str("");
            filename << SOURCE_DIR << "/res/integration/test_semi_euler_smooth_off"
                     << steps << ".vtk";

            solver.integrationStep();

            save_particles_to_vtk(filename.str(), solver.getSystem().getPositions(),
                                  solver.getSystem().getDensities());

            std::cout << "Results saved to " << filename.str() << std::endl;            
        }
    }
}

TEST_CASE("SolverRun", "[simulation]") {
    std::cout << "Testing time constant simulation.." << std::endl;

    const double particleDiameter = 0.1;

    // Sample Particles in a Box
    FluidSystem particles = Emitter().sampleFluidBox(Eigen::Vector3d(0.0, 0.25, 0),
                                                     Eigen::Vector3d(1.0, 1.25, 1.0),
                                                     particleDiameter);

    SolverSPH solver(particles);
    solver.setSnapShotAfterMS(40);
    solver.setParamStiffness(0.0);
    solver.setFluidViscosity(0.0);
    solver.enableGravity(false);
    solver.enableSmoothing(false);

    SECTION("SimpleSolverI") {      
        solver.run("solver_test_I", 12000);
    }

    solver.setParamStiffness(1000);
    
    SECTION("SimpleSolverII") {        
        solver.run("solver_test_II", 12000);
    }

    solver.addBoundary(System::Emitter().sampleBoundaryPlane(
        Eigen::Vector3d(0.0, 0.0, 0.0),
        Eigen::Vector3d(1.0, 0.0, 0.0),
        Eigen::Vector3d(0.0, 0.0, 1.0),
        Eigen::Vector3d(1.0, 0.0, 1.0),
        particleDiameter));
    
    solver.enableGravity(true);
    
    SECTION("SimpleSolverIII") {       
        solver.run("solver_test_III", 12000);
    }

    solver.setFluidViscosity(0.1);
    
    SECTION("SimpleSolverIV") {
        solver.run("solver_test_IV", 6000);
    }

    SECTION("SimpleSolverV") {
        solver.setBoundaryViscosity(0, 1000);
        solver.run("solver_test_V", 6000);
    }

    solver.enableSmoothing(true);
    
    SECTION("SimpleSolverVI") {
        solver.enableGravity(true);
        solver.run("solver_test_VI", 6000);
    }
}

TEST_CASE("ComplexRun", "[complex]") {
    SECTION("DAMNBREAK") {
        std::cout << "Testing complex simulation.." << std::endl;

        const double particleDiameter = 0.05;

        // Sample Particles in a Box
        FluidSystem particles = System::Emitter().sampleFluidBox(
            Eigen::Vector3d(0, 0.0, 0),
            Eigen::Vector3d(1.0, 1.0, 1.0),
            particleDiameter);

        SolverSPH solver(particles);
        solver.setSnapShotAfterMS(40);
        solver.setParamStiffness(1000.0);
        solver.setFluidViscosity(0.02);
        solver.enableGravity(true);
        solver.enableSmoothing(true);

        solver.addBoundary(System::Emitter().sampleBoundaryHollowBox(
            Eigen::Vector3d(-0.05, -0.05, -0.05),
            Eigen::Vector3d(2.5, 2.5, 2.5),
            particleDiameter));
        solver.setBoundaryViscosity(0, 0.02);
                                                      
        solver.run("complex_simulation", 6000);
    }
}
