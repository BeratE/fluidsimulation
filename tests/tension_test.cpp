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
        particles.setGamma(0.01);
        
        Solver *solver;
        std::stringstream filename;
        filename << "tension";
        
        SECTION("SPH") {
            solver = new SolverSPH(particles);
            ((SolverSPH*)solver)->setParamStiffness(3000);

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
        solver->setMaxTimeStepSeconds(0.002);
        solver->enableGravity(false);
        solver->enableSmoothing(true);
        solver->enableAdhesion(false);

        solver->run(filename.str(), 2000);       
    }        
}


TEST_CASE("Adhesion", "[ahesion]") {
    SECTION("Adhesion") {
        const double particleDiameter = 0.08;

        FluidSystem particles = Emitter().sampleFluidBox(Eigen::Vector3d(-.75, 1.5, -.75),
                                                         Eigen::Vector3d( .75, 3.5,  .75),
                                                         particleDiameter);

        particles.setGamma(0.05);
        
        std::stringstream filepath;
        filepath << SOURCE_DIR << "/res/" << "icosphere.obj";
        BoundarySystem icosphere =
            System::Emitter().sampleBoundaryMesh(filepath.str(), particleDiameter);

        icosphere.setBeta(20);
        
        Solver *solver;
        std::stringstream filename;
        filename << "adhesion";
        
        SECTION("SPH") {
            solver = new SolverSPH(particles);
            ((SolverSPH*)solver)->setParamStiffness(3000);

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
        solver->setMaxTimeStepSeconds(0.002);
        solver->enableGravity(true);
        solver->enableSmoothing(true);
        solver->enableTension(false);

        solver->run(filename.str(), 2000);       
    }
}



TEST_CASE("Tension on funnel", "[tension_funnel]") {
    //omp_set_dynamic(0);
    //omp_set_num_threads(1);
    SECTION("Funnel") {
        std::cout << "Funnel with tension and adhesion" << std::endl;

        const double particleDiameter = 0.05;
        const double particelDiameterBoundary = particleDiameter / 0.5;

        FluidSystem particles = Emitter().sampleFluidBox(Eigen::Vector3d(0, 0, 0),
                                                         Eigen::Vector3d(1, 1, 1),
                                                         particleDiameter);
        particles.setGamma(0.2);

        BoundarySystem plane1 = Emitter().sampleBoundaryPlane(Eigen::Vector3d(-1, 0, -1),
                                                              Eigen::Vector3d(0.3, -1, 0.3),
                                                              Eigen::Vector3d(2, 0, -1),
                                                              Eigen::Vector3d(0.7, -1, 0.3),
                                                              particelDiameterBoundary);
        //plane1.setViscosity(0.02);

        BoundarySystem plane2 = Emitter().sampleBoundaryPlane(Eigen::Vector3d(2, 0, -1),
                                                              Eigen::Vector3d(0.7, -1, 0.3),
                                                              Eigen::Vector3d(2, 0, 2),
                                                              Eigen::Vector3d(0.7, -1, 0.7),
                                                              particelDiameterBoundary);
        //plane2.setViscosity(0.02);

        BoundarySystem plane3 = Emitter().sampleBoundaryPlane(Eigen::Vector3d(2, 0, 2),
                                                              Eigen::Vector3d(0.7, -1, 0.7),
                                                              Eigen::Vector3d(-1, 0, 2),
                                                              Eigen::Vector3d(0.3, -1, 0.7),
                                                              particelDiameterBoundary);
        //plane3.setViscosity(0.02);

        BoundarySystem plane4 = Emitter().sampleBoundaryPlane(Eigen::Vector3d(-1, 0, 2),
                                                              Eigen::Vector3d(0.3, -1, 0.7),
                                                              Eigen::Vector3d(-1, 0, -1),
                                                              Eigen::Vector3d(0.3, -1, 0.3),
                                                              particelDiameterBoundary);
        //plane4.setViscosity(0.02);

        BoundarySystem box = Emitter().sampleBoundaryHollowBox(Eigen::Vector3d(-1, -1.5, 2),
                                                               Eigen::Vector3d(2, 1.2, -1),
                                                               particelDiameterBoundary);
        box.setViscosity(0.08);

        SECTION("SPH") {
            SolverSPH solver(particles);
            solver.addBoundary(plane1);
            solver.addBoundary(plane2);
            solver.addBoundary(plane3);
            solver.addBoundary(plane4);
            solver.addBoundary(box);

            solver.setSnapShotAfterMS(40);
            solver.setParamStiffness(3000);
            solver.setFluidViscosity(0.1);
            solver.enableAdhesion(false);
            solver.enableGravity(true);
            solver.enableSmoothing(true);
            solver.enableTension(false);

            SECTION("OFF") {
                solver.run("sph_funnel_no_tension_no_adhesion", 3000);
            }

            SECTION("BETA250") {
                solver.enableAdhesion(true);
                solver.enableTension(true);
                solver.run("sph_funnel_1_beta", 3000);
            }
        }

        SECTION("PBF") {
            SolverPBF solver(particles);
            solver.addBoundary(plane1);
            solver.addBoundary(plane2);
            solver.addBoundary(plane3);
            solver.addBoundary(plane4);
            solver.addBoundary(box);

            solver.setParamSmoothing(0.25);
            solver.setSnapShotAfterMS(40);
            solver.setFluidViscosity(0.0);
            solver.enableAdhesion(false);
            solver.enableGravity(true);
            solver.enableSmoothing(true);
            solver.enableTension(false);
            solver.setNumIterations(3);

            SECTION("OFF") {
                solver.run("pbf_funnel_no_tension_no_adhesion", 3000);
            }

            SECTION("BETA250") {
                solver.enableAdhesion(true);
                solver.enableTension(true);
                solver.run("pbf_funnel_1_beta", 3000);
            }
        }
    }
}

TEST_CASE("Surface Tension Complex Scene", "[tension_complex]") {
    std::cout << "Complex scene with tension and adhesion" << std::endl;

    const double particleDiameter = 0.05;
    const double particelDiameterBoundary = particleDiameter / 0.5;

    FluidSystem particles = Emitter().sampleFluidBox( Eigen::Vector3d(0.1, 1.1, 1.9),//Eigen::Vector3d(0.05, 1.05, 1.95),
                                                      Eigen::Vector3d(1.9, 1.9, 0.1),
                                                      particleDiameter);
    particles.setGamma(0.2);

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
        SolverSPH solver(particles);
        solver.addBoundary(bigBoundary);
        solver.addBoundary(plane1);
        solver.addBoundary(plane2);
        solver.addBoundary(plane3);
        solver.addBoundary(plane4);
        solver.addBoundary(armadillo);

        solver.setSnapShotAfterMS(40);
        solver.setParamStiffness(3000);
        solver.setFluidViscosity(0.1);
        solver.enableAdhesion(false);
        solver.enableGravity(true);
        solver.enableSmoothing(true);
        solver.enableTension(false);

        SECTION("SPH_tension_complex_no_tension") {
            solver.run("sph_tension_complex_no_tension", 3000);
        }

        SECTION("SPH_tension_complex_no_tension") {
            solver.enableAdhesion(true);
            solver.enableTension(true);
            solver.run("sph_tension_complex", 3000);
        }
    }


    SECTION("PBF") {
        SolverPBF solver(particles);
        solver.addBoundary(bigBoundary);
        solver.addBoundary(plane1);
        solver.addBoundary(plane2);
        solver.addBoundary(plane3);
        solver.addBoundary(plane4);
        solver.addBoundary(armadillo);

        solver.setParamSmoothing(0.25);
        solver.setSnapShotAfterMS(40);
        solver.setFluidViscosity(0.0);
        solver.enableAdhesion(false);
        solver.enableGravity(true);
        solver.enableSmoothing(true);
        solver.enableTension(false);
        solver.setNumIterations(3);

        SECTION("PBF_tension_complex") {
            solver.run("pbf_tension_complex_no_tension", 3000);
        }

        SECTION("SPH_tension_complex_no_tension") {
            solver.enableAdhesion(true);
            solver.enableTension(true);
            solver.run("pbf_tension_complex", 3000);
        }
    }
}
