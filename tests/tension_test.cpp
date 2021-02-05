#include "catch.hpp"
#include "config.h"
#include "vtk_writer.h"
#include "learnSPH/solver_sph.h"
#include "learnSPH/solver_pbf.h"
#include "learnSPH/system/fluidsystem.h"
#include "learnSPH/system/emitter.h"
#include <omp.h>

using namespace learnSPH;
using namespace learnSPH::System;

TEST_CASE("Simple surface tension, no gravity cube") {
   SECTION("Tension") {
       // omp_set_dynamic(0);
       // omp_set_num_threads(1);        
       
       const double particleDiameter = 0.08;

       FluidSystem particles = Emitter().sampleFluidBox(Eigen::Vector3d(0, 0, 0),
                                                        Eigen::Vector3d(1, 1, 1),
                                                        particleDiameter);
       particles.setGamma(0.1);
       
       Solver *solver;
       std::stringstream filename;
       filename << SOURCE_DIR << "/res/simulation/" << "tension_two";
       
       SECTION("SPH") {
           solver = new SolverSPH(particles);
           ((SolverSPH*)solver)->setParamStiffness(3000);
           solver->setMaxTimeStepSeconds(0.001);
           
           filename << "_sph";
           SECTION("ON") {
               solver->enableTension(true);
           }
           
           SECTION("OFF") {
               solver->enableTension(false);
               filename << "_off";
           }
       }                
       SECTION("PBF") {
           solver = new SolverPBF(particles);
           ((SolverPBF*)solver)->setNumIterations(3);
           solver->setMaxTimeStepSeconds(0.002);
           filename << "pbf";
           
           SECTION("ON") {
               solver->enableTension(true);
           }
           SECTION("OFF") {
               solver->enableTension(false);
               filename << "_off";
           }
       }

       solver->setSnapShotAfterMS(33);
       solver->setFluidViscosity(0.001);
       
       solver->enableGravity(false);
       solver->enableSmoothing(true);
       solver->enableAdhesion(false);

       double startTime = omp_get_wtime();
       
       solver->run(filename.str(), 3000);
       double endTime = omp_get_wtime();

       std::cout << "Runtime: " << endTime-startTime << std::endl;
   }        
}


TEST_CASE("Adhesion", "[adhesion]") {
    SECTION("Adhesion") {
        double startTime = omp_get_wtime();
        
        const double particleDiameter = 0.1;

        Eigen::Vector3d origin(0, 2.0, 0);
        auto sphereSDF = [origin](Eigen::Vector3d x)
            { return (x - origin).norm() - 0.5; };

        FluidSystem particles = Emitter().sampleFluidSDF(sphereSDF,
                                                         Eigen::Vector3d(-0.5,  1.5, -0.5),
                                                         Eigen::Vector3d(   1,  1,  1),
                                                         particleDiameter);
        
        // FluidSystem particles = Emitter().sampleFluidBox(Eigen::Vector3d(-.75, 1.5, -.75),
        //                                                  Eigen::Vector3d( .75, 3.5,  .75),
        //                                                  particleDiameter);

        particles.setGamma(0.1);
        
        std::stringstream filepath;
        filepath << SOURCE_DIR << "/res/" << "icosphere.obj";
        BoundarySystem icosphere =
            System::Emitter().sampleBoundaryMesh(filepath.str(), particleDiameter);

        icosphere.setBeta(2);
        
        Solver *solver;
        std::stringstream filename;
        filename << SOURCE_DIR << "/res/simulation/" << "adhesion_two";
        
        SECTION("SPH") {
            solver = new SolverSPH(particles);
            ((SolverSPH*)solver)->setParamStiffness(3000);
            solver->setMaxTimeStepSeconds(0.001);

            filename << "_sph";
            SECTION("ON") {
                solver->enableTension(true);
                solver->enableAdhesion(true);
            }
            
            SECTION("OFF") {
                solver->enableAdhesion(false);
                filename << "_off";
            }
        }                
        SECTION("PBF") {
            solver = new SolverPBF(particles);
            ((SolverPBF*)solver)->setNumIterations(3);
            solver->setMaxTimeStepSeconds(0.002);
            filename << "pbf";
            
            SECTION("ON") {
                solver->enableTension(true);
                solver->enableAdhesion(true);
            }
            SECTION("OFF") {
                solver->enableAdhesion(false);
                filename << "_off";
            }
        }

        solver->addBoundary(icosphere);
        solver->setSnapShotAfterMS(33);
        solver->setFluidViscosity(0.01);
        solver->setBoundaryViscosity(0, 0.01);
        solver->enableGravity(true);
        solver->enableSmoothing(true);

        solver->run(filename.str(), 2000);

        double endTime = omp_get_wtime();

        std::cout << "Runtime: " << endTime-startTime << std::endl;
    }
}

//////////////////////////////////////////////////////////////////////////
/// Bitte nochmal auf die Tests schauen
//////////////////////////////////////////////////////////////////////////

TEST_CASE("Surface Tension Complex Scene", "[tension_complex]") {
    std::cout << "Complex scene with tension and adhesion" << std::endl;

    const double particleDiameter = 0.05;
    const double particelDiameterBoundary = particleDiameter / 0.5;

    FluidSystem particles = Emitter().sampleFluidBox(Eigen::Vector3d(0.1, 1.1, 1.9),
        Eigen::Vector3d(1.9, 1.9, 0.1),
        particleDiameter);

    BoundarySystem bigBoundary = Emitter().sampleBoundaryHollowBox(Eigen::Vector3d(0.0, 0.0, 0.0),
        Eigen::Vector3d(4.0, 2.0, 2.0),
        particelDiameterBoundary);
    bigBoundary.setViscosity(0.02);

    BoundarySystem plane1 = Emitter().sampleBoundaryPlane(Eigen::Vector3d(0.0, 1.0, 2.0),
        Eigen::Vector3d(2.0, 1.0, 2.0),
        Eigen::Vector3d(0.0, 1.0, 0.0),
        Eigen::Vector3d(2.0, 1.0, 0.0),
        particelDiameterBoundary);
    plane1.setViscosity(0.02);

    BoundarySystem plane2 = Emitter().sampleBoundaryPlane(Eigen::Vector3d(2.0, 2.0, 2.0),
        Eigen::Vector3d(2.0, 1.5, 2.0),
        Eigen::Vector3d(2.0, 2.0, 0.0),
        Eigen::Vector3d(2.0, 1.5, 0.0),
        particelDiameterBoundary);
    plane2.setViscosity(0.02);

    BoundarySystem plane3 = Emitter().sampleBoundaryPlane(Eigen::Vector3d(2.0, 1.5, 2.0),
        Eigen::Vector3d(2.0, 1.0, 2.0),
        Eigen::Vector3d(2.0, 1.5, 1.25),
        Eigen::Vector3d(2.0, 1.0, 1.25),
        particelDiameterBoundary);
    plane3.setViscosity(0.02);

    BoundarySystem plane4 = Emitter().sampleBoundaryPlane(Eigen::Vector3d(2.0, 1.5, 0.75),
        Eigen::Vector3d(2.0, 1.0, 0.75),
        Eigen::Vector3d(2.0, 1.5, 0.0),
        Eigen::Vector3d(2.0, 1.0, 0.0),
        particelDiameterBoundary);
    plane4.setViscosity(0.02);

    std::stringstream filepath;
    filepath << SOURCE_DIR << "/res/" << "armadillo_low.obj";
    BoundarySystem armadillo = Emitter().sampleBoundaryMesh(filepath.str(), particelDiameterBoundary);
    armadillo.setViscosity(0.08);

    Eigen::Matrix4d transform = Eigen::Matrix4d::Identity();
    transform(0, 3) = 2.5;
    transform(1, 3) = 0.1;
    transform(2, 3) = 1.0;
    armadillo.transform(transform);

    SECTION("SPH") {
        std::stringstream filename;
        filename << SOURCE_DIR << "/res/simulation/" << "sph_complex";
        particles.setGamma(0.5); // Gamma etwas niedriger, siehe Beispielvideo von Jose
        SolverSPH solver(particles);

        solver.setSnapShotAfterMS(40);
        solver.setParamStiffness(3000);
        solver.setFluidViscosity(0.07);
        solver.enableAdhesion(false);
        solver.enableGravity(true);
        solver.enableSmoothing(true);
        solver.enableTension(false);

        SECTION("SPH_tension_complex_no_tension") {
            filename << "_no_tension";
        }

        SECTION("SPH_tension_complex_no_tension") {
            solver.enableAdhesion(true);
            solver.enableTension(true);
            armadillo.setBeta(5);
            filename << "_5_beta";
        }

        solver.addBoundary(bigBoundary);
        solver.addBoundary(plane1);
        solver.addBoundary(plane2);
        solver.addBoundary(plane3);
        solver.addBoundary(plane4);
        solver.addBoundary(armadillo);
        solver.run(filename.str(), 3000);
    }


    SECTION("PBF") {
        std::stringstream filename;
        filename << SOURCE_DIR << "/res/simulation/" << "pbf_complex";
        particles.setGamma(0.1);
        SolverPBF solver(particles);

        solver.setParamSmoothing(0.25);
        solver.setSnapShotAfterMS(40);
        solver.setFluidViscosity(0.01);
        solver.enableAdhesion(false);
        solver.enableGravity(true);
        solver.enableSmoothing(true);
        solver.enableTension(false);
        solver.setNumIterations(3);

        SECTION("SPH_tension_complex_no_tension") {
            filename << "_no_tension";
        }

        SECTION("SPH_tension_complex_no_tension") {
            solver.enableAdhesion(true);
            solver.enableTension(true);
            armadillo.setBeta(5);
            filename << "_5_beta";
        }

        solver.addBoundary(bigBoundary);
        solver.addBoundary(plane1);
        solver.addBoundary(plane2);
        solver.addBoundary(plane3);
        solver.addBoundary(plane4);
        solver.addBoundary(armadillo);
        solver.run(filename.str(), 3000);
    }
}