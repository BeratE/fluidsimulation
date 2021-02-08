#include "catch.hpp"
#include "config.h"
#include "vtk_writer.h"
#include "learnSPH/kernel.h"
#include "learnSPH/solver.h"
#include "learnSPH/solver_sph.h"
#include "learnSPH/solver_pbf.h"
#include "learnSPH/system/emitter.h"
#include "learnSPH/system/particlesystem.h"
#include "learnSPH/system/boundarysystem.h"
#include "surface/surface.h"
#include <omp.h>

using namespace learnSPH;
using namespace learnSPH::System;
using namespace learnSPH::Surface;

// TEST_CASE("Armadillo Splash","")
// {    
//    SECTION("Armadillo_Splash") {
       
//        const double particleDiameter = 0.025;
//        const double boundaryDiameter = particleDiameter * 1.6;

//        // Sample Particles in a Box
//        FluidSystem particles = System::Emitter().sampleFluidBox(
//            Eigen::Vector3d(-0.5, 1.2, -0.5),
//            Eigen::Vector3d( 0.5, 1.8,  0.5),
//            particleDiameter);

//        std::stringstream filepath;
//        filepath << SOURCE_DIR << "/res/" << "armadillo_low.obj";
       
//        BoundarySystem armadillo = System::Emitter().sampleBoundaryMesh(
//            filepath.str(),
//            boundaryDiameter);
       
//        BoundarySystem box = System::Emitter().sampleBoundaryHollowBox(
//            Eigen::Vector3d(-1.05, -0.10, -1.05),
//            Eigen::Vector3d( 1.05,  2.00,  1.05),
//            particleDiameter);        

//        std::stringstream filename;
//        std::stringstream outputfile;
//        Solver *solver;
       
//        SECTION("SPH_I") {
//             solver = new SolverSPH(particles);
//             filename << "_sph_I";
            
//             ((SolverSPH*)solver)->setParamStiffness(1000);            
//             solver->setMaxTimeStepSeconds(0.0001);
//             solver->setParamSmoothing(0.25);

//             solver->setFluidViscosity(0.005);
//             solver->setFluidTension(0.1);

//             solver->addBoundary(icosphere);
//             solver->setBoundaryAdhesion(0, 500); // icosphere
//             solver->setBoundaryViscosity(0, 0.02);

//             solver->addBoundary(funnel);
//             solver->setBoundaryAdhesion(1, 0.01); // funnel
//             solver->setBoundaryViscosity(1, 0.005);
//        }
//        SECTION("PBF_V") {
//             solver = new SolverPBF(particles);
//             filename << "_pbf_V";

//             ((SolverPBF*)solver)->setNumIterations(3);
//             solver->setMaxTimeStepSeconds(0.014);
//             solver->setParamSmoothing(0.1);

//             solver->setFluidViscosity(0.0001);
//             solver->setFluidTension(0.15);

//             solver->addBoundary(icosphere);
//             solver->setBoundaryAdhesion(0, 100); // icosphere
//             solver->setBoundaryViscosity(0, 0.0008);

//             solver->addBoundary(funnel);
//             solver->setBoundaryAdhesion(1, 0.5); // funnel
//             solver->setBoundaryViscosity(1, 0.0005);
//         }              

//        solver->setParamDrag(0.1);
//        solver->setSnapShotAfterMS(1000.0 / 40);
//        solver->enableGravity(true);
//        solver->enableSmoothing(true);
//        solver->enableTension(true);
//        solver->enableAdhesion(true);
       
//        solver->addBoundary(box);
//        solver->addBoundary(armadillo);

//        solver->setBoundaryViscosity(0, 0.02);
//        solver->setBoundaryViscosity(1, 0.08);
//        outputfile << "Armadillo Viscosity: " << 0.08 << std::endl;
//        outputfile << "Box Viscosity: " << 0.02 << std::endl;

//        const double start_t = omp_get_wtime();
       
//        solver->run(filename.str(), 6000);
       
//        const double end_t = omp_get_wtime();
//        const double delta_t = end_t - start_t;
       

//        outputfile << "Num Particles: " << particles.getSize() << std::endl;
//        outputfile << "Runtime: " << delta_t << std::endl;

//        std::string tmp = filename.str();
//        filename.str(std::string());
//        filename << SOURCE_DIR << "/res/simulation/" << tmp << ".out";
       
//        std::ofstream outfile;
//        outfile.open(filename.str());
//        outfile << outputfile.rdbuf();
//        outfile.close();

//        delete solver;        
//    }
// }

TEST_CASE("SmallDam") {
    SECTION("SmallDam") {
        const double particleDiameter = 0.025;
        const double boundaryDiameter = particleDiameter;
        FluidSystem particles = System::Emitter().sampleFluidBox(
            Eigen::Vector3d(0, 0, 0),
            Eigen::Vector3d(0.5, 1.45, 1),
            particleDiameter);

        BoundarySystem box = System::Emitter().sampleBoundaryHollowBox(
            Eigen::Vector3d(-0.05, -0.05, -0.05), Eigen::Vector3d(1.00, 1.5, 1.05),
            boundaryDiameter);

        Solver* solver;

        std::stringstream filename;
        filename << SOURCE_DIR << "/res/simulation2/" << "smallDam";

        SECTION("SPH") {
            filename << "_SPH";
            solver = new SolverSPH(particles);
            solver->addBoundary(box);

            solver->enableGravity(true);
            solver->enableAdhesion(true);
            solver->enableDrag(true);
            solver->enableSmoothing(true);
            solver->enableTension(true);

            SECTION("I") {
                solver->setParamDrag(0.025);
                filename << "_I";
                ((SolverSPH*)solver)->setParamStiffness(3000);
                solver->setMaxTimeStepSeconds(0.001);
                solver->setParamSmoothing(0.5);

                solver->setFluidViscosity(0.001);
                solver->enableTension(false);

                solver->setBoundaryViscosity(0, 0.001); // box
                solver->enableAdhesion(false);
            }

            SECTION("II") {
                solver->setParamDrag(0.025);
                filename << "_II";
                ((SolverSPH*)solver)->setParamStiffness(3000);
                solver->setMaxTimeStepSeconds(0.001);
                solver->setParamSmoothing(0.5);

                solver->setFluidViscosity(0.025);
                solver->enableTension(false);

                solver->setBoundaryViscosity(0, 0.001); // box
                solver->enableAdhesion(false);
            }

            SECTION("III") {
                solver->setParamDrag(0.025);
                filename << "_III";
                ((SolverSPH*)solver)->setParamStiffness(3000);
                solver->setMaxTimeStepSeconds(0.001);
                solver->setParamSmoothing(0.5);

                solver->setFluidViscosity(0.001);
                solver->enableTension(false);

                solver->setBoundaryViscosity(0, 0.025); // box
                solver->enableAdhesion(false);
            }

            SECTION("IV") {
                solver->setParamDrag(0.025);
                filename << "_IV";
                ((SolverSPH*)solver)->setParamStiffness(3000);
                solver->setMaxTimeStepSeconds(0.001);
                solver->setParamSmoothing(0.5);

                solver->setFluidViscosity(0.001);
                solver->setFluidTension(0.25);

                solver->setBoundaryViscosity(0, 0.001); // box
                solver->setBoundaryAdhesion(0, 10);
            }

            SECTION("V") {
                solver->setParamDrag(0.05);
                filename << "_V";
                ((SolverSPH*)solver)->setParamStiffness(3000);
                solver->setMaxTimeStepSeconds(0.001);
                solver->setParamSmoothing(0.5);

                solver->setFluidViscosity(0.001);
                solver->enableTension(false);

                solver->setBoundaryViscosity(0, 0.001); // box
                solver->enableAdhesion(false);
            }

            SECTION("VI") {
                solver->setParamDrag(0.05);
                filename << "_VI";
                ((SolverSPH*)solver)->setParamStiffness(3000);
                solver->setMaxTimeStepSeconds(0.001);
                solver->setParamSmoothing(0.5);

                solver->setFluidViscosity(0.025);
                solver->enableTension(false);

                solver->setBoundaryViscosity(0, 0.001); // box
                solver->enableAdhesion(false);
            }

            SECTION("VII") {
                solver->setParamDrag(0.05);
                filename << "_VII";
                ((SolverSPH*)solver)->setParamStiffness(3000);
                solver->setMaxTimeStepSeconds(0.001);
                solver->setParamSmoothing(0.5);

                solver->setFluidViscosity(0.001);
                solver->enableTension(false);

                solver->setBoundaryViscosity(0, 0.025); // box
                solver->enableAdhesion(false);
            }

            SECTION("VIII") {
                solver->setParamDrag(0.05);
                filename << "_VIII";
                ((SolverSPH*)solver)->setParamStiffness(3000);
                solver->setMaxTimeStepSeconds(0.001);
                solver->setParamSmoothing(0.5);

                solver->setFluidViscosity(0.001);
                solver->setFluidTension(0.25);

                solver->setBoundaryViscosity(0, 0.001); // box
                solver->setBoundaryAdhesion(0, 10);
            }
        }

        SECTION("PBF") {
            filename << "_PBF";
            solver = new SolverPBF(particles);
            solver->addBoundary(box);

            solver->enableGravity(true);
            solver->enableAdhesion(true);
            solver->enableDrag(true);
            solver->enableSmoothing(true);
            solver->enableTension(true);

            SECTION("I") {
                solver->setParamDrag(0.025);
                filename << "_I";
                ((SolverPBF*)solver)->setNumIterations(7);
                solver->setMaxTimeStepSeconds(0.004);
                solver->setParamSmoothing(0.1);

                solver->setFluidViscosity(0.001);
                solver->enableTension(false);

                solver->setBoundaryViscosity(0, 0.001); // box
                solver->enableAdhesion(false);
            }

            SECTION("II") {
                solver->setParamDrag(0.025);
                filename << "_II";
                ((SolverPBF*)solver)->setNumIterations(7);
                solver->setMaxTimeStepSeconds(0.004);
                solver->setParamSmoothing(0.1);

                solver->setFluidViscosity(0.025);
                solver->enableTension(false);

                solver->setBoundaryViscosity(0, 0.001); // box
                solver->enableAdhesion(false);
            }

            SECTION("III") {
                solver->setParamDrag(0.025);
                filename << "_III";
                ((SolverPBF*)solver)->setNumIterations(7);
                solver->setMaxTimeStepSeconds(0.004);
                solver->setParamSmoothing(0.1);

                solver->setFluidViscosity(0.001);
                solver->enableTension(false);

                solver->setBoundaryViscosity(0, 0.025); // box
                solver->enableAdhesion(false);
            }

            SECTION("IV") {
                solver->setParamDrag(0.025);
                filename << "_IV";
                ((SolverPBF*)solver)->setNumIterations(7);
                solver->setMaxTimeStepSeconds(0.004);
                solver->setParamSmoothing(0.1);

                solver->setFluidViscosity(0.001);
                solver->setFluidTension(0.25);

                solver->setBoundaryViscosity(0, 0.001); // box
                solver->setBoundaryAdhesion(0, 10);
            }

            SECTION("V") {
                solver->setParamDrag(0.05);
                filename << "_V";
                ((SolverPBF*)solver)->setNumIterations(7);
                solver->setMaxTimeStepSeconds(0.004);
                solver->setParamSmoothing(0.1);

                solver->setFluidViscosity(0.001);
                solver->enableTension(false);

                solver->setBoundaryViscosity(0, 0.001); // box
                solver->enableAdhesion(false);
            }

            SECTION("VI") {
                solver->setParamDrag(0.05);
                filename << "_VI";
                ((SolverPBF*)solver)->setNumIterations(7);
                solver->setMaxTimeStepSeconds(0.004);
                solver->setParamSmoothing(0.1);

                solver->setFluidViscosity(0.025);
                solver->enableTension(false);

                solver->setBoundaryViscosity(0, 0.001); // box
                solver->enableAdhesion(false);
            }

            SECTION("VII") {
                solver->setParamDrag(0.05);
                filename << "_VII";
                ((SolverPBF*)solver)->setNumIterations(7);
                solver->setMaxTimeStepSeconds(0.004);
                solver->setParamSmoothing(0.1);

                solver->setFluidViscosity(0.001);
                solver->enableTension(false);

                solver->setBoundaryViscosity(0, 0.004); // box
                solver->enableAdhesion(false);
            }

            SECTION("VIII") {
                solver->setParamDrag(0.05);
                filename << "_VIII";
                ((SolverPBF*)solver)->setNumIterations(7);
                solver->setMaxTimeStepSeconds(0.004);
                solver->setParamSmoothing(0.1);

                solver->setFluidViscosity(0.001);
                solver->setFluidTension(0.25);

                solver->setBoundaryViscosity(0, 0.001); // box
                solver->setBoundaryAdhesion(0, 10);
            }
        }

        solver->setSnapShotAfterMS(1000.0 / 40);

        const double startTime = omp_get_wtime();

        solver->run(filename.str(), 4000);

        double endTime = omp_get_wtime();

        std::cout << "Runtime: " << endTime - startTime << std::endl;

        (filename.str(), *solver, endTime - startTime);

        delete solver;
    }
    
}

TEST_CASE("ArmadilloBreak")
{
   SECTION("ArmadilloArmy") {       
       const double particleDiameter = 0.03;
       const double boundaryDiameter = particleDiameter;
       
       // Sample Particles in a Box
       FluidSystem particles = System::Emitter().sampleFluidBox(
           Eigen::Vector3d(-1.0, 0.0, -1.2), Eigen::Vector3d(1.0, 2.1, -0.2),
           particleDiameter);
       std::cout << "Num Particles: " << particles.getSize() << std::endl;

       BoundarySystem box = System::Emitter().sampleBoundaryHollowBox(
           Eigen::Vector3d(-1.05, -0.05, -1.25), Eigen::Vector3d(1.05, 2.15, 3.5),
           boundaryDiameter * 1.1);

       std::stringstream filepath;
       filepath << SOURCE_DIR << "/res/" << "armadillo_low.obj";
       BoundarySystem armadillo =
           System::Emitter().sampleBoundaryMesh(
               filepath.str(), boundaryDiameter);
            
       Solver *solver;

       std::stringstream filename;
       filename << SOURCE_DIR << "/res/armadillo/" << "armadilloArmy";

       SECTION("SPH") {
           solver = new SolverSPH(particles);
           filename << "_sph";

           solver->addBoundary(box);
           
           Eigen::Matrix4d transform;       
           transform = Eigen::Matrix4d::Identity();
           transform(2, 3) = 0.5;
           armadillo.transform(transform);        
           solver->addBoundary(armadillo);           

           transform(2, 3) = 1.0;
           transform(0, 3) = 0.5;
           armadillo.transform(transform);
           solver->addBoundary(armadillo);
       
           transform(0, 3) = -1.0;
           armadillo.transform(transform);
           solver->addBoundary(armadillo);

           SECTION("V") {
               filename << "_V";
               ((SolverSPH *)solver)->setParamStiffness(2000);
               solver->setMaxTimeStepSeconds(0.002);
               solver->setParamSmoothing(0.25);

               solver->setFluidViscosity(0.004);
               solver->setFluidTension(0.5);

               solver->setBoundaryViscosity(0, 0.004); // box
               solver->setBoundaryAdhesion(0, 0.8);

               for (int i = 1; i < 4; i++) { // armadillos
                   solver->setBoundaryViscosity(i, 0.004 * i);
                   solver->setBoundaryAdhesion(i, 0.5 * i);
               }
           }

           SECTION("VI") {
               filename << "_VI";
               ((SolverSPH*)solver)->setParamStiffness(2000);
               solver->setMaxTimeStepSeconds(0.001);
               solver->setParamSmoothing(0.5);

               solver->setFluidViscosity(0.004);
               solver->setFluidTension(0.5);

               solver->setBoundaryViscosity(0, 0.004); // box
               solver->setBoundaryAdhesion(0, 0.8);

               for (int i = 1; i < 4; i++) { // armadillos
                   solver->setBoundaryViscosity(i, 0.004 * i);
                   solver->setBoundaryAdhesion(i, 0.5 * i);
               }
           }
           
           SECTION("VII") {
               filename << "_VII";
               ((SolverSPH*)solver)->setParamStiffness(2000);
               solver->setMaxTimeStepSeconds(0.001);
               solver->setParamSmoothing(0.25);

               solver->setFluidViscosity(0.004);
               solver->setFluidTension(0.25);

               solver->setBoundaryViscosity(0, 0.004); // box
               solver->setBoundaryAdhesion(0, 0.8);

               for (int i = 1; i < 4; i++) { // armadillos
                   solver->setBoundaryViscosity(i, 0.004 * i);
                   solver->setBoundaryAdhesion(i, 0.5 * i);
               }
           }

           SECTION("VIII") {
               filename << "_VIII";
               ((SolverSPH*)solver)->setParamStiffness(2000);
               solver->setMaxTimeStepSeconds(0.001);
               solver->setParamSmoothing(0.5);

               solver->setFluidViscosity(0.004);
               solver->setFluidTension(0.25);

               solver->setBoundaryViscosity(0, 0.004); // box
               solver->setBoundaryAdhesion(0, 0.8);

               for (int i = 1; i < 4; i++) { // armadillos
                   solver->setBoundaryViscosity(i, 0.004 * i);
                   solver->setBoundaryAdhesion(i, 0.5 * i);
               }
           }
           SECTION("II") {
               filename << "_II";
               ((SolverSPH *)solver)->setParamStiffness(3000);
               solver->setMaxTimeStepSeconds(0.0005);
               solver->setParamSmoothing(0.50);

               solver->setFluidViscosity(0.01);
               solver->setFluidTension(0.25);

               solver->setBoundaryViscosity(0, 0.008); // box
               solver->setBoundaryAdhesion(0, 1.4);

               for (int i = 1; i < 4; i++) { // armadillos
                   solver->setBoundaryViscosity(i, 0.04 * i);
                   solver->setBoundaryAdhesion(i, 2.0 * i);
               }
           }
       }

       SECTION("PBF") {
           solver = new SolverPBF(particles);
           filename << "_pbf";

           solver->addBoundary(box);
           
           Eigen::Matrix4d transform;       
           transform = Eigen::Matrix4d::Identity();
           transform(2, 3) = 0.5;
           armadillo.transform(transform);        
           solver->addBoundary(armadillo);           

           transform(2, 3) = 1.0;
           transform(0, 3) = 0.5;
           armadillo.transform(transform);
           solver->addBoundary(armadillo);
       
           transform(0, 3) = -1.0;
           armadillo.transform(transform);
           solver->addBoundary(armadillo);

           SECTION("III") {
               filename << "_III";
               ((SolverPBF *)solver)->setNumIterations(3);
               solver->setMaxTimeStepSeconds(0.014);
               solver->setParamSmoothing(0.1);

               solver->setFluidViscosity(0.0001);
               solver->setFluidTension(0.1);

               solver->setBoundaryViscosity(0, 0.0004); // box
               solver->setBoundaryAdhesion(0, 0.2);
               
               for (int i = 1; i < 4; i++) {  // armadillos
                   solver->setBoundaryViscosity(i, 0.0004 * i);
                   solver->setBoundaryAdhesion(i, 0.5 * i);
               }                      
           }
           SECTION("IV") {
               filename << "_IV";
               ((SolverPBF*)solver)->setNumIterations(3);
               solver->setMaxTimeStepSeconds(0.002);
               solver->setParamSmoothing(0.1);

               solver->setFluidViscosity(0.0001);
               solver->setFluidTension(0.1);

               solver->setBoundaryViscosity(0, 0.001); // box
               solver->setBoundaryAdhesion(0, 1.0);
               
               for (int i = 1; i < 4; i++) {  // armadillos
                   solver->setBoundaryViscosity(i, 0.004 * i);
                   solver->setBoundaryAdhesion(i, 6.0 * i);
               }                      
           }
       }                     

       solver->setParamDrag(0.05);
       solver->setSnapShotAfterMS(1000.0 / 40);
       solver->enableGravity(true);
       solver->enableSmoothing(true);
       solver->enableTension(true);
       solver->enableAdhesion(true);              

       const double startTime = omp_get_wtime();
       
       solver->run(filename.str(), 4000);

       double endTime = omp_get_wtime();

       std::cout << "Runtime: " << endTime - startTime << std::endl;

       outputParams(filename.str(), *solver, endTime - startTime);

       delete solver;    
   }
}
