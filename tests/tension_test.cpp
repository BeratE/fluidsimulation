#include "catch.hpp"
#include "config.h"
#include "vtk_writer.h"
#include "learnSPH/solver_sph.h"
#include "learnSPH/solver_pbf.h"
#include "learnSPH/system/fluidsystem.h"
#include "learnSPH/system/emitter.h"
#include <omp.h>
#include "util.h"

using namespace learnSPH;
using namespace learnSPH::System;


TEST_CASE("Simple surface tension, no gravity cube") {
  SECTION("Tension") {
      //omp_set_dynamic(0);
      //omp_set_num_threads(1);        
      
      const double particleDiameter = 0.025;

      FluidSystem particles = Emitter().sampleFluidBox(Eigen::Vector3d(0, 0, 0),
                                                       Eigen::Vector3d(1, 1, 1),
                                                       particleDiameter);
      
      Solver *solver;
      std::stringstream filename;
      filename << SOURCE_DIR << "/res/big_simulations/" << "floating_ball";
      
      SECTION("SPH") {
          solver = new SolverSPH(particles);
          filename << "_gamma_1_smoothing_0,5_eta_1,2_sph";
          
          ((SolverSPH*)solver)->setParamStiffness(3000);          
          solver->setMaxTimeStepSeconds(0.001);

          solver->setFluidTension(1.0); 
          solver->setFluidViscosity(0.01);          
      }                
      SECTION("PBF") {
          solver = new SolverPBF(particles);
          filename << "_gamma_1_smoothing_0,25_eta_1,2_pbf";
          
          ((SolverPBF*)solver)->setNumIterations(3);          
          solver->setMaxTimeStepSeconds(0.002);
          solver->setParamSmoothing(0.25);
          solver->setFluidTension(1.0);
          solver->setFluidViscosity(0.0008);
      }

      solver->setSnapShotAfterMS(1000.0/60.0);
      
      solver->enableGravity(false);
      solver->enableSmoothing(true);
      solver->enableAdhesion(false);
      solver->enableTension(true);

      double startTime = omp_get_wtime();

      solver->run(filename.str(), 2500);

      double endTime = omp_get_wtime();

      std::cout << "Runtime: " << endTime-startTime << std::endl;

      outputParams(filename.str(), *solver, endTime-startTime);

      delete solver;
  }        
}


TEST_CASE("Adhesion", "[adhesion]") {
    SECTION("Adhesion") {
        double startTime = omp_get_wtime();

        const double particleDiameter = 0.05;
        const double boundaryDiameter = particleDiameter * 1.4;

        FluidSystem particles = Emitter().sampleFluidBox(
            Eigen::Vector3d(-0.8, 4.0, -.4),
            Eigen::Vector3d(-1.8, 4.8, .4),
            particleDiameter);

        std::stringstream filepath;
        filepath << SOURCE_DIR << "/res/" << "icosphere.obj";
        BoundarySystem icosphere =
            System::Emitter().sampleBoundaryMesh(filepath.str(), boundaryDiameter);


        filepath.str(std::string());
        filepath << SOURCE_DIR << "/res/" << "funnel.obj";
        BoundarySystem funnel =
            System::Emitter().sampleBoundaryMesh(filepath.str(), boundaryDiameter);   

        Solver* solver;
        std::stringstream filename;
        filename << SOURCE_DIR << "/res/simulation2/" << "waterOnBall";
        
        SECTION("SPH") {
            solver = new SolverSPH(particles);
            filename << "_sph";
            
            ((SolverSPH*)solver)->setParamStiffness(3000);            
            solver->setMaxTimeStepSeconds(0.0005);
            
            solver->setFluidViscosity(0.01);
            solver->setFluidTension(0.2);

            solver->addBoundary(icosphere);
            solver->setBoundaryAdhesion(0, 2); // icosphere
            solver->setBoundaryViscosity(0, 0.02);

            solver->addBoundary(funnel);
            solver->setBoundaryAdhesion(1, 0.5); // funnel
            solver->setBoundaryViscosity(1, 0.005);
        }
        SECTION("PBF") {
            solver = new SolverPBF(particles);
            filename << "_pbf";

            ((SolverPBF*)solver)->setNumIterations(5);
            solver->setMaxTimeStepSeconds(0.0005);
            
            solver->setFluidViscosity(0.0001);
            solver->setFluidTension(0.05);

            solver->addBoundary(icosphere);
            solver->setBoundaryAdhesion(0, 1); // icosphere
            solver->setBoundaryViscosity(0, 0.008);

            solver->addBoundary(funnel);
            solver->setBoundaryAdhesion(1, 0.4); // funnel
            solver->setBoundaryViscosity(1, 0.0005);                       
        }        
            
        solver->setSnapShotAfterMS(1000.0 / 40);
        solver->enableGravity(true);
        solver->enableSmoothing(true);
        solver->enableTension(false);
        solver->enableAdhesion(false);

        solver->run(filename.str(), 1500);

        double endTime = omp_get_wtime();

        std::cout << "Runtime: " << endTime - startTime << std::endl;

        outputParams(filename.str(), *solver, endTime - startTime);

        delete solver;
    }
}


TEST_CASE("Water droplet") {
    SECTION("WaterDroplet") {

        auto sceneSDF = [](Eigen::Vector3d x)
            {
                Eigen::Vector3d boxOrigin(0, 0.3, 0);
                Eigen::Vector3d boxBounds(4, 0.6, 4);
                auto boxSDF = [boxOrigin, boxBounds](Eigen::Vector3d x) {
                    x = x - boxOrigin;
                    Eigen::Vector3d q =
                        Eigen::Vector3d(x.array().abs()) - boxBounds / 2;
                    return Eigen::Vector3d(q.array().max(0.0)).norm() +
                        std::min(std::max(q(0), std::max(q(1), q(2))), 0.0);
                };

                Eigen::Vector3d sphereOrigin(0, 2.0, 0);
                double sphereRadius = 0.6;
                auto sphereSDF = [sphereOrigin, sphereRadius](Eigen::Vector3d x) {
                    return (x - sphereOrigin).norm() - sphereRadius;
                };

                return std::min(sphereSDF(x), boxSDF(x));
            };

        double particleRadius = 0.1;
        FluidSystem particles = Emitter().sampleFluidSDF(sceneSDF,
                                     Eigen::Vector3d(-2, 0, -2),
                                     Eigen::Vector3d( 4, 4,  4),
                                     particleRadius);
        particles.setViscosity(0.01);
        particles.setGamma(0.15);


        BoundarySystem box = Emitter().sampleBoundaryHollowBox(
            Eigen::Vector3d(-2.2, -0.2, -2.2),
            Eigen::Vector3d( 2.2,  4.2,  2.2),
            particleRadius);
        box.setViscosity(0.02);
        box.setBeta(0.5);
        

       Solver *solver;
       std::stringstream filename;
       filename << SOURCE_DIR << "/res/simulation2/" << "waterdroplet";
       
       SECTION("SPH") {
           filename << "_sph";
           solver = new SolverSPH(particles);
           ((SolverSPH*)solver)->setParamStiffness(3000);
           solver->setMaxTimeStepSeconds(0.001);
           solver->setFluidViscosity(0.01);
       }                
       SECTION("PBF") {
           filename << "_pbf";
           solver = new SolverPBF(particles);
           ((SolverPBF*)solver)->setNumIterations(3);
           solver->setMaxTimeStepSeconds(0.002);
           solver->setFluidViscosity(0.001);
       }

       solver->addBoundary(box);

       solver->setSnapShotAfterMS(1000.0/60.0);
       solver->enableGravity(true);
       solver->enableSmoothing(true);
       solver->enableAdhesion(true);
       solver->enableTension(true);

       double startTime = omp_get_wtime();              

       solver->run(filename.str(), 5000);

       double endTime = omp_get_wtime();

       std::cout << "Runtime: " << endTime-startTime << std::endl;

       delete solver;
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
        filename << SOURCE_DIR << "/res/simulation2/" << "sph_complex";
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
        filename << SOURCE_DIR << "/res/simulation2/" << "pbf_complex";
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
