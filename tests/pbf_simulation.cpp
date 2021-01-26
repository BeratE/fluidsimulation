//#include "catch.hpp"
//#include "config.h"
//#include "vtk_writer.h"
//#include "learnSPH/kernel.h"
//#include "learnSPH/solver_pbf.h"
//#include "learnSPH/system/emitter.h"
//#include "learnSPH/system/particlesystem.h"
//
//using namespace learnSPH;
//using namespace learnSPH::System;
//
//TEST_CASE("PBF Solver Run", "[pbf]") {
//    SECTION("PBFSolver") {
//        const double particleDiameter = 0.1;
//
//        // Sample Particles in a Box
//        FluidSystem particles = Emitter().sampleFluidBox(Eigen::Vector3d(0.0, 0.25, 0),
//                                                         Eigen::Vector3d(1.0, 1.25, 1.0),
//                                                         particleDiameter);
//
//        SolverPBF solver(particles);
//        solver.setSnapShotAfterMS(40);
//        solver.setFluidViscosity(0.0);
//        solver.enableSmoothing(false);
//        solver.enableGravity(false);
//   
//        SECTION("PBFSolverI") {      
//            solver.run("pbf_solver_test_I", 6000);
//        }
//
//
//    
//        SECTION("PBFSolverII") {      
//            solver.run("pbf_solver_test_II", 6000);
//        }
//
//        solver.enableGravity(true);
//    
//        solver.addBoundary(System::Emitter().sampleBoundaryPlane(
//                               Eigen::Vector3d(0.0, 0.0, 0.0),
//                               Eigen::Vector3d(1.0, 0.0, 0.0),
//                               Eigen::Vector3d(0.0, 0.0, 1.0),
//                               Eigen::Vector3d(1.0, 0.0, 1.0),
//                               particleDiameter));
//    
//        SECTION("PBFSolverIII") {       
//            solver.run("pbf_solver_test_III", 6000);
//        }
//
//        solver.setFluidViscosity(0.1);
//    
//        SECTION("PBFSolverIV") {
//            solver.run("pbf_solver_test_IV", 6000);
//        }
//
//        SECTION("PBFSolverV") {
//            solver.setBoundaryViscosity(0, 1000);
//            solver.run("pbf_solver_test_V", 6000);
//        }
//
//        solver.enableSmoothing(true);
//
//        SECTION("PBFSolverVI") {
//            solver.run("pbf_solver_test_VI", 6000);
//        }
//    }
//}
//
//TEST_CASE("ComplexPBF", "")
//{
//    SECTION("PBFDambreak") {
//        const double particleDiameter = 0.05;
//
//        // Sample Particles in a Box
//        FluidSystem particles = System::Emitter().sampleFluidBox(
//            Eigen::Vector3d(0, 0.0, 0),
//            Eigen::Vector3d(1.0, 1.0, 1.0),
//            particleDiameter);
//
//        SolverPBF solver(particles);
//        solver.setSnapShotAfterMS(40);
//        solver.setFluidViscosity(0.01);
//        solver.enableGravity(true);
//        solver.enableSmoothing(true);
//
//        solver.addBoundary(System::Emitter().sampleBoundaryHollowBox(
//                               Eigen::Vector3d(-0.05, -0.05, -0.05),
//                               Eigen::Vector3d(2.5, 2.5, 2.5),
//                               particleDiameter));
//        solver.setBoundaryViscosity(0, 0.02);
//                                                      
//        solver.run("pbf_complex_simulation", 6000);
//    }
//}
