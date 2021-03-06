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

//TEST_CASE("Simple surface tension, no gravity cube") {
//  SECTION("Tension") {
//      // omp_set_dynamic(0);
//      // omp_set_num_threads(1);        
//      
//      const double particleDiameter = 0.035;
//
//      FluidSystem particles = Emitter().sampleFluidBox(Eigen::Vector3d(0, 0, 0),
//                                                       Eigen::Vector3d(1, 1, 1),
//                                                       particleDiameter);
//      
//      Solver *solver;
//      std::stringstream filename;
//      filename << SOURCE_DIR << "/res/simulation2/" << "floating_ball";
//      
//      SECTION("SPH") {
//          filename << "_sph";
//          solver = new SolverSPH(particles);
//          
//          ((SolverSPH*)solver)->setParamStiffness(3000);
//          solver->setMaxTimeStepSeconds(0.002);
//
//          solver->setFluidViscosity(0.01);
//          solver->setFluidTension(0.15);
//      }                
//      SECTION("PBF") {
//          filename << "_pbf";
//          solver = new SolverPBF(particles);
//          
//          ((SolverPBF*)solver)->setNumIterations(3);
//          solver->setMaxTimeStepSeconds(0.01);
//          
//          solver->setFluidViscosity(0.001);
//          solver->setFluidTension(0.15);
//      }
//
//      solver->setSnapShotAfterMS(1000.0/30.0);
//      
//      solver->enableGravity(false);
//      solver->enableSmoothing(false);
//      solver->enableAdhesion(false);
//      solver->enableTension(true);
//
//      double startTime = omp_get_wtime();
//
//      solver->run(filename.str(), 3000);
//
//      double endTime = omp_get_wtime();
//
//      std::cout << "Runtime: " << endTime-startTime << std::endl;
//
//      outputParams(filename.str(), *solver, endTime-startTime);
//
//      delete solver;
//  }        
//}


TEST_CASE("Adhesion", "[adhesion]") {
    SECTION("Adhesion") {
        
        const double particleDiameter = 0.025;
        const double boundaryDiameter = particleDiameter;

        FluidSystem particles = Emitter().sampleFluidBox(
            Eigen::Vector3d(-0.8, 4.0, -.4),
            Eigen::Vector3d(-1.8, 4.8, .4),
            particleDiameter);
        std::cout << "Num Particles: " << particles.getSize() << std::endl;

        std::stringstream filepath;
        filepath << SOURCE_DIR << "/res/" << "icosphere.obj";
        BoundarySystem icosphere =
            System::Emitter().sampleBoundaryMesh(filepath.str(), boundaryDiameter);

        Eigen::Matrix4d transform = Eigen::Matrix4d::Identity();
        transform(0, 3) = -0.3 - 1.2886002659484357;
        transform(1, 3) = 0.3 + 0.9638809794597053;
        transform(2, 3) = 0.0 - 0.044312775519577456;
        icosphere.transform(transform);

        /*filepath.str(std::string());
        filepath << SOURCE_DIR << "/res/" << "funnel.obj";
        BoundarySystem funnel =
            System::Emitter().sampleBoundaryMesh(filepath.str(), boundaryDiameter);   */

        Solver* solver;
        std::stringstream filename;
        filename << SOURCE_DIR << "/res/waterOnBall/" << "waterOnBall";
       
        
        SECTION("SPH") {
            solver = new SolverSPH(particles);
            

            ((SolverSPH*)solver)->setParamStiffness(1000);
            solver->setMaxTimeStepSeconds(0.0005);
            solver->setParamSmoothing(0.25);

            solver->setFluidViscosity(0.005);

            solver->addBoundary(icosphere);
            solver->setBoundaryAdhesion(0, 250);
            
            SECTION("SPH_XIV") {
                filename << "_sph_XIV";
                solver->enableTension(true);
                solver->setFluidTension(0.01);
            }

            SECTION("SPH_XV") {
                filename << "_sph_XV";
                solver->enableTension(true);
                solver->setFluidTension(0.025);
            }

            SECTION("SPH_XVI") {
                filename << "_sph_XVI";
                solver->enableTension(false);
            }
            solver->setBoundaryViscosity(0, 0.02);

            // solver->addBoundary(funnel);
            // solver->setBoundaryAdhesion(1, 0.01); // funnel
            // solver->setBoundaryViscosity(1, 0.005);
        }


        SECTION("PBF") {

            solver = new SolverPBF(particles);

            ((SolverPBF*)solver)->setNumIterations(3);
            solver->setMaxTimeStepSeconds(0.002);
            solver->setParamSmoothing(0.1);
            
            solver->setFluidViscosity(0.001);

            solver->addBoundary(icosphere);
            solver->setBoundaryAdhesion(0, 250); // icosphere

            SECTION("PBF_XV") {
                filename << "_pbf_XV";
                solver->enableTension(true);
                solver->setFluidTension(0.01);
            }
            SECTION("PBF_XVI") {
                filename << "_pbf_XVI";
                solver->enableTension(true);
                solver->setFluidTension(0.025);
            }
            SECTION("PBF_XVII") {
                filename << "_pbf_XVII";
                solver->enableTension(false);
            }
            solver->setBoundaryViscosity(0, 0.008);

            // solver->addBoundary(funnel);
            // solver->setBoundaryAdhesion(1, 0.5); // funnel
            // solver->setBoundaryViscosity(1, 0.0005);
        }
        
        solver->setParamDrag(0.1);
        solver->setSnapShotAfterMS(1000.0 / 40);
        solver->enableGravity(true);
        solver->enableSmoothing(true);
        solver->enableAdhesion(true);
        
        double startTime = omp_get_wtime();
        
        solver->run(filename.str(), 4000);

        double endTime = omp_get_wtime();
        outputParams(filename.str(), *solver, endTime - startTime);

        std::cout << "Runtime: " << endTime - startTime << std::endl;


        delete solver;
    }
}


TEST_CASE("Water droplet") {
    SECTION("WaterDroplet") {

        auto sceneSDF = [](Eigen::Vector3d x)
            {
                Eigen::Vector3d boxOrigin(0.0, 0.25, 0.0);
                Eigen::Vector3d boxBounds(2.5, 0.50, 2.5);
                auto boxSDF = [boxOrigin, boxBounds](Eigen::Vector3d x) {
                    x = x - boxOrigin;
                    Eigen::Vector3d q =
                        Eigen::Vector3d(x.array().abs()) - boxBounds / 2;
                    return Eigen::Vector3d(q.array().max(0.0)).norm() +
                        std::min(std::max(q(0), std::max(q(1), q(2))), 0.0);
                };

                double sphereRadius = 0.4;
                Eigen::Vector3d sphereOrigin(0, 1.6, 0);
                auto sphereSDF = [sphereOrigin, sphereRadius](Eigen::Vector3d x) {
                    return (x - sphereOrigin).norm() - sphereRadius;
                };

                return std::min(sphereSDF(x), boxSDF(x));
            };

        double particleDiameter = 0.042;
	    double boundaryDiameter = particleDiameter * 1.6;
        FluidSystem particles = Emitter().sampleFluidSDF(sceneSDF,
                                     Eigen::Vector3d(-2, 0, -2),
                                     Eigen::Vector3d( 4, 4,  4),
                                     particleDiameter);

        std::cout << "Num Particles: " << particles.getSize() << std::endl;

        BoundarySystem box = Emitter().sampleBoundaryHollowBox(
            Eigen::Vector3d(-1.32, -0.08, -1.32),
            Eigen::Vector3d( 1.32,  2.20,  1.32),
            boundaryDiameter);
        

       Solver *solver;
       std::stringstream filename;
       filename << SOURCE_DIR << "/res/simulation2/" << "waterdroplet";
       
       SECTION("SPH_I") {
           filename << "_sph_I";
           solver = new SolverSPH(particles);

           ((SolverSPH*)solver)->setParamStiffness(3000);
           solver->setMaxTimeStepSeconds(0.002);

	   solver->setFluidTension(0.20);
           solver->setFluidViscosity(0.001);

	   solver->addBoundary(box);
	   solver->setBoundaryAdhesion(0, 0.8); // box
	   solver->setBoundaryViscosity(0, 0.002);
	   
       }
       SECTION("SPH_II") {
           filename << "_sph_II";
           solver = new SolverSPH(particles);           

           ((SolverSPH*)solver)->setParamStiffness(5000);
           solver->setMaxTimeStepSeconds(0.002);
           solver->setParamSmoothing(0.25);

	   solver->setFluidTension(0.50);
           solver->setFluidViscosity(0.001);

	   solver->addBoundary(box);
	   solver->setBoundaryAdhesion(0, 0.2); // box
	   solver->setBoundaryViscosity(0, 0.001);	   
       }
       SECTION("SPH_III") {
           filename << "_sph_III";
           solver = new SolverSPH(particles);           

           ((SolverSPH*)solver)->setParamStiffness(3000);
           solver->setMaxTimeStepSeconds(0.001);
           solver->setParamSmoothing(0.25);

	   solver->setFluidTension(1.0);
           solver->setFluidViscosity(0.01);

	   solver->addBoundary(box);
	   solver->setBoundaryAdhesion(0, 1.0); // box
	   solver->setBoundaryViscosity(0, 0.001);	   
       }
       SECTION("SPH_IV") {
           filename << "_sph_IIV";
           solver = new SolverSPH(particles);           

           ((SolverSPH*)solver)->setParamStiffness(4000);
           solver->setMaxTimeStepSeconds(0.0005);
           solver->setParamSmoothing(0.25);

	   solver->setFluidTension(1.0);
           solver->setFluidViscosity(0.01);

	   solver->addBoundary(box);
	   solver->setBoundaryAdhesion(0, 0.8); // box
	   solver->setBoundaryViscosity(0, 0.001);	   
       }
       SECTION("PBF_I") {
           filename << "_pbf_I";
           solver = new SolverPBF(particles);

           ((SolverPBF*)solver)->setNumIterations(4);
           solver->setMaxTimeStepSeconds(0.012);
           solver->setParamSmoothing(0.1);

	   solver->setFluidTension(0.05);
           solver->setFluidViscosity(0.0001);

	   solver->addBoundary(box);
	   solver->setBoundaryAdhesion(0, 0.2); // box
	   solver->setBoundaryViscosity(0, 0.0001);
       }
       SECTION("PBF_II") {
           filename << "_pbf_II";
           solver = new SolverPBF(particles);

           ((SolverPBF*)solver)->setNumIterations(5);
           solver->setMaxTimeStepSeconds(0.028);
           solver->setParamSmoothing(0.1);

	   solver->setFluidTension(0.15);
           solver->setFluidViscosity(0.0001);

	   solver->addBoundary(box);
	   solver->setBoundaryAdhesion(0, 0.2); // box
	   solver->setBoundaryViscosity(0, 0.0001);
       }
       SECTION("PBF_III") {
           filename << "_pbf_III";
           solver = new SolverPBF(particles);

           ((SolverPBF*)solver)->setNumIterations(6);
           solver->setMaxTimeStepSeconds(0.016);
           solver->setParamSmoothing(0.1);

	   solver->setFluidTension(0.15);
           solver->setFluidViscosity(0.001);

	   solver->addBoundary(box);
	   solver->setBoundaryAdhesion(0, 0.2); // box
	   solver->setBoundaryViscosity(0, 0.001);
       }
       SECTION("PBF_IV") {
           filename << "_pbf_IV";
           solver = new SolverPBF(particles);

           ((SolverPBF*)solver)->setNumIterations(5);
           solver->setMaxTimeStepSeconds(0.004);
           solver->setParamSmoothing(0.1);

	   solver->setFluidTension(1.0);
           solver->setFluidViscosity(0.01);

	   solver->addBoundary(box);
	   solver->setBoundaryAdhesion(0, 1.0); // box
	   solver->setBoundaryViscosity(0, 0.01);
       }

	     
       solver->setParamDrag(0.01);

       solver->setSnapShotAfterMS(1000.0/60.0);
       solver->enableGravity(true);
       solver->enableSmoothing(true);
       solver->enableAdhesion(true);
       solver->enableTension(true);

       double startTime = omp_get_wtime();              

       solver->run(filename.str(), 4000);

       double endTime = omp_get_wtime();

       std::cout << "Runtime: " << endTime-startTime << std::endl;
       
       outputParams(filename.str(), *solver, endTime - startTime);
       
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
