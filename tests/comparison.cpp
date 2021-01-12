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

TEST_CASE("Comparison","")
{    
    SECTION("Armadillo_Splash") {
        double ratioSmoothingLengthSamplingStep = 2.0;
        const double particleDiameter = 0.05;

        // Sample Particles in a Box
        FluidSystem particles = System::Emitter().sampleFluidBox(
            Eigen::Vector3d(-0.5, 1.2, -0.5),
            Eigen::Vector3d( 0.5, 1.8,  0.5),
            particleDiameter);

        std::stringstream filepath;
        filepath << SOURCE_DIR << "/res/" << "armadillo_low.obj";
        
        BoundarySystem armadillo = System::Emitter().sampleBoundaryMesh(
            filepath.str(),
            particleDiameter);
        
        BoundarySystem box = System::Emitter().sampleBoundaryHollowBox(
            Eigen::Vector3d(-1.05, -0.10, -1.05),
            Eigen::Vector3d( 1.05,  2.00,  1.05),
            particleDiameter);        


        std::stringstream filename;
        std::stringstream outputfile;
        Solver *solver;
        
        SECTION("SPH") {
            double stiffness;
            double viscosity;
            bool smoothing;
            
            SECTION("ArmadilloSPH_I") {                
                filename << "armadillo_sph_I";
                stiffness = 1000;
                viscosity = 0.01;
                smoothing = true;                
            }

            SECTION("ArmadilloSPH_II") {                
                filename << "armadillo_sph_II";
                stiffness = 3000;
                viscosity = 0.1;
                smoothing = true;                
            }

            SECTION("ArmadilloSPH_III") {                
                filename << "armadillo_sph_III";
                stiffness = 5000;
                viscosity = 0.4;
                smoothing = true;                
            }
            
            solver = new SolverSPH(particles);
            ((SolverSPH *)solver)->setParamStiffness(stiffness);
            solver->setFluidViscosity(viscosity);
            solver->enableSmoothing(smoothing);

            outputfile << filename.str() << std::endl
                       << "Stiffness: " << stiffness << std::endl
                       << "Smoothing: " << smoothing << std::endl
                       << "Fluid Viscosity: " << viscosity << std::endl;
        }

        SECTION("PBF") {
            size_t iterations;
            double viscosity;
            bool smoothing;
            
            SECTION("ArmadilloPBF_I") {
                filename << "armadillo_pbf_I";
                iterations = 3;
                viscosity = 0.1;
                smoothing = true;
            }

            SECTION("ArmadilloPBF_II") {
                filename << "armadillo_pbf_II";
                iterations = 3;
                viscosity = 0.0;
                smoothing = true;
            }

            SECTION("ArmadilloPBF_III") {
                filename << "armadillo_pbf_III";
                iterations = 5;
                viscosity = 0.2;
                smoothing = false;
            }

            solver = new SolverPBF(particles);
            ((SolverPBF *)solver)->setNumIterations(iterations);
            solver->setFluidViscosity(viscosity);
            solver->enableSmoothing(smoothing);

            outputfile << filename.str() << std::endl
                       << "Iterations: " << iterations << std::endl
                       << "Smoothing:" << smoothing << std::endl
                       << "Fluid Viscosity: " << viscosity << std::endl;
        }
        
        solver->setSnapShotAfterMS(30);
        solver->enableGravity(true);
        solver->addBoundary(box);
        solver->addBoundary(armadillo);

        solver->setBoundaryViscosity(0, 0.02);
        solver->setBoundaryViscosity(1, 0.08);
        outputfile << "Armadillo Viscosity: " << 0.08 << std::endl;
        outputfile << "Box Viscosity: " << 0.02 << std::endl;

        const double start_t = omp_get_wtime();
        
        solver->run(filename.str(), 6000);
        
        const double end_t = omp_get_wtime();
        const double delta_t = end_t - start_t;
        

        outputfile << "Num Particles: " << particles.getSize() << std::endl;
        outputfile << "Runtime: " << delta_t << std::endl;

        std::string tmp = filename.str();
        filename.str(std::string());
        filename << SOURCE_DIR << "/res/simulation/" << tmp << ".out";
        
        std::ofstream outfile;
        outfile.open(filename.str());
        outfile << outputfile.rdbuf();
        outfile.close();

        delete solver;        
    }
}

TEST_CASE("BigCompare")
{
    SECTION("ArmadilloArmy") {
        double ratioSmoothingLengthSamplingStep = 2.0;
        const double particleDiameter = 0.05;

        // Sample Particles in a Box
        FluidSystem particles = System::Emitter().sampleFluidBox(
            Eigen::Vector3d(-1.0, 0.2, -1.0), Eigen::Vector3d(1.0, 1.8, -0.1),
            particleDiameter);
             
        std::stringstream filename;
        std::stringstream outputfile;
        Solver *solver;

        SECTION("SPH") {
            double stiffness;
            double viscosity;
            bool smoothing;            

            SECTION("ArmadilloSPH_I") {                
                filename << "2complex_armadillo_sph_I";
                stiffness = 3000;
                viscosity = 0.04;
                smoothing = true;                
            }
            
            solver = new SolverSPH(particles);
            ((SolverSPH *)solver)->setParamStiffness(stiffness);
            solver->setFluidViscosity(viscosity);
            solver->enableSmoothing(smoothing);

            outputfile << filename.str() << std::endl
                       << "Stiffness: " << stiffness << std::endl
                       << "Smoothing: " << smoothing << std::endl
                       << "Fluid Viscosity: " << viscosity << std::endl;
        }
        SECTION("PBF") {
            size_t iterations;
            double viscosity;
            bool smoothing;
            
            SECTION("ArmadilloPBF_I") {
                filename << "complex_armadillo_pbf_I";
                iterations = 4;
                viscosity = 0.1;
                smoothing = false;
            }

            solver = new SolverPBF(particles);
            ((SolverPBF *)solver)->setNumIterations(iterations);
            solver->setFluidViscosity(viscosity);
            solver->enableSmoothing(smoothing);

            outputfile << filename.str() << std::endl
                       << "Iterations: " << iterations << std::endl
                       << "Smoothing:" << smoothing << std::endl
                       << "Fluid Viscosity: " << viscosity << std::endl;
        }

        solver->setSnapShotAfterMS(30);
        solver->enableGravity(true);

        BoundarySystem box = System::Emitter().sampleBoundaryHollowBox(
            Eigen::Vector3d(-1.05, -0.10, -1.05), Eigen::Vector3d(1.05, 2.00, 3.05),
            particleDiameter);
        
        solver->addBoundary(box);
        solver->setBoundaryViscosity(0, 0.02);

        std::stringstream filepath;
        filepath << SOURCE_DIR << "/res/" << "armadillo_low.obj";
        BoundarySystem armadillo =
            System::Emitter().sampleBoundaryMesh(filepath.str(), particleDiameter);

        Eigen::Matrix4d transform;
        
        transform = Eigen::Matrix4d::Identity();
        transform(2, 3) = 0.5;
        armadillo.transform(transform);        
        solver->addBoundary(armadillo);
        solver->setBoundaryViscosity(1, 0.04);

        transform(2, 3) = 1.0;
        transform(0, 3) = 0.5;
        armadillo.transform(transform);
        solver->addBoundary(armadillo);
        solver->setBoundaryViscosity(2, 0.04);
        
        transform(0, 3) = -1.0;
        armadillo.transform(transform);
        solver->addBoundary(armadillo);
        solver->setBoundaryViscosity(3, 0.04);
        
        
        outputfile << "Armadillo Viscosity: " << 0.04 << std::endl;
        outputfile << "Box Viscosity: " << 0.02 << std::endl;

        const double start_t = omp_get_wtime();
        
        solver->run(filename.str(), 12000);
        
        const double end_t = omp_get_wtime();
        const double delta_t = end_t - start_t;
        

        outputfile << "Num Particles: " << particles.getSize() << std::endl;
        outputfile << "Runtime: " << delta_t << std::endl;

        std::string tmp = filename.str();
        filename.str(std::string());
        filename << SOURCE_DIR << "/res/simulation/" << tmp << ".out";
        
        std::ofstream outfile;
        outfile.open(filename.str());
        outfile << outputfile.rdbuf();
        outfile.close();

        delete solver;        
    }
}
