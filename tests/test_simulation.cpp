#include "catch.hpp"
#include "particlesystem.h"
#include "emitter.h"
#include "kernel.h"
#include "solver.h"
#include "util/vtk_writer.h"
#include "util/config.h"

using namespace learnSPH;

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
    FluidSystem particles = Emitter().sampleFluidBox(Eigen::Vector3d(0, 0.0, 0),
                                                     Eigen::Vector3d(1.0, 1.0, 1.0),
                                                     particleDiameter);

    SolverSPH solver(particles);
    solver.setSnapShotAfterMS(40);
    solver.setParameterStiffness(0.0);
    solver.setParameterViscosity(0.0);
    solver.enableGravity(false);
    solver.enableSmoothing(false);

    SECTION("SimpleSolverI") {      
        solver.run("solver_test_I", 6000);
    }
    SECTION("SimpleSolverII") {
        solver.setParameterStiffness(1000);
        solver.run("solver_test_II", 6000);
    }
    SECTION("SimpleSolverIII") {
        solver.enableGravity(true);
        solver.setParameterStiffness(1000);
        solver.addBoundary(Emitter().sampleBoundaryPlane(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                         Eigen::Vector3d(1.0, 0.0, 0.0),
                                                         Eigen::Vector3d(0.0, 0.0, 1.0),
                                                         Eigen::Vector3d(1.0, 0.0, 1.0),
                                                         particleDiameter));        
        solver.run("solver_test_III", 6000);
    }
    
    SECTION("SimpleSolverIV") {
        solver.enableGravity(true);
        solver.setParameterStiffness(1000);
        BoundarySystem plane = Emitter().sampleBoundaryPlane(Eigen::Vector3d(0.0, 0.0, 0.0),
            Eigen::Vector3d(1.0, 0.0, 0.0),
            Eigen::Vector3d(0.0, 0.0, 1.0),
            Eigen::Vector3d(1.0, 0.0, 1.0),
            particleDiameter);
        solver.addBoundary(plane);
        solver.setParameterViscosity(0.1);
        solver.run("solver_test_IV", 6000);
    }

    SECTION("SimpleSolverV") {
        solver.enableGravity(true);
        solver.setParameterStiffness(1000);
        BoundarySystem plane = Emitter().sampleBoundaryPlane(Eigen::Vector3d(0.0, 0.0, 0.0),
            Eigen::Vector3d(1.0, 0.0, 0.0),
            Eigen::Vector3d(0.0, 0.0, 1.0),
            Eigen::Vector3d(1.0, 0.0, 1.0),
            particleDiameter);
        plane.setViscosity(1000);
        solver.addBoundary(plane);
        solver.setParameterViscosity(0.1);


        solver.run("solver_test_V", 6000);
    }

    SECTION("SimpleSolverVI") {
        solver.enableGravity(true);
        solver.setParameterStiffness(1000);
        solver.enableSmoothing(true);
        BoundarySystem plane = Emitter().sampleBoundaryPlane(Eigen::Vector3d(0.0, 0.0, 0.0),
            Eigen::Vector3d(1.0, 0.0, 0.0),
            Eigen::Vector3d(0.0, 0.0, 1.0),
            Eigen::Vector3d(1.0, 0.0, 1.0),
            particleDiameter);
        solver.addBoundary(plane);

        solver.run("solver_test_VI", 6000);
    }
}

TEST_CASE("ComplexRun", "[complex]") {
    SECTION("DAMNBREAK") {
        std::cout << "Testing complex simulation.." << std::endl;

        const double particleDiameter = 0.05;

        // Sample Particles in a Box
        FluidSystem particles = Emitter().sampleFluidBox(
                                                         Eigen::Vector3d(0, 0.0, 0), Eigen::Vector3d(1.0, 1.0, 1.0),
                                                         particleDiameter);

        SolverSPH solver(particles);
        solver.setSnapShotAfterMS(40);
        solver.setParameterStiffness(1000.0);
        solver.setParameterViscosity(0.02);
        solver.enableGravity(true);
        solver.enableSmoothing(true);

        BoundarySystem box = Emitter().sampleBoundaryHollowBox(
                                                               Eigen::Vector3d(-0.05, -0.05, -0.05), Eigen::Vector3d(2.5, 2.5, 2.5),
                                                               particleDiameter);
        box.setViscosity(0.02);
        solver.addBoundary(box);
        solver.run("complex_simulation", 6000);
    }
}