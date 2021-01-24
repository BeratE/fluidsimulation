#include "catch.hpp"
#include <tuple>
#include <string>
#include <iostream>
#include <Eigen/Dense>
#include "vtk_writer.h"
#include "config.h"
#include "surface/surface.h"
#include <omp.h>
#include "learnSPH/system/emitter.h"
#include "learnSPH/system/particlesystem.h"
#include "learnSPH/kernel.h"
#include "learnSPH/solver_sph.h"
#include <CompactNSearch/CompactNSearch.h>

using namespace learnSPH;
using namespace learnSPH::Surface;
using namespace learnSPH::System;



TEST_CASE("Construction", "") {
    //printCudaVersion();
    SECTION("Construction") {
        
        const double start_t = omp_get_wtime();

        const Eigen::Vector3i numVerts(30, 30, 30);
        const Eigen::Vector3d gridSize(1.0, 1.0, 1.0);
    
        std::vector<Eigen::Vector3d> gridVerts;
        std::vector<double> gridSDF;
        bool run = false;
    
        SECTION("SphereSDF") {
            const double radius = 0.2;
            const Eigen::Vector3d origin(0.5, 0.5, 0.5);
            auto sdf = [origin,radius](Eigen::Vector3d x)
                { return (x-origin).norm() - radius; };
        
            discretizeSDF(gridSize, numVerts, sdf,
                          &gridSDF, &gridVerts);
            run = true;
        }
        SECTION("TorusSDF") {
            const double r = 0.1;
            const double R = 0.3;
            const Eigen::Vector3d origin(0.5, 0.5, 0.5);
            auto sdf = [r, R, origin](Eigen::Vector3d x)
                {
                    x = x - origin;
                    return r*r-pow(sqrt(x(0)*x(0)+x(1)*x(1))-R, 2)-x(2)*x(2);
                };
        
            discretizeSDF(gridSize, numVerts, sdf,
                          &gridSDF, &gridVerts);
            run = true;
        }
        if (run) {
            std::vector<Eigen::Vector3d> vertices;
            std::vector<std::array<int, 3>> triangles;
            marchCubes(numVerts, gridSDF, gridVerts, vertices, triangles);

            std::stringstream filename;
            filename << SOURCE_DIR << "/res/surface/simple_surface.vtk";

            save_mesh_to_vtk(filename.str(), vertices, triangles);
            // save_particles_to_vtk(filename.str(), verts);
        }
        const double end_t = omp_get_wtime();
        const double delta_t = end_t - start_t;

        std::cout << "Runtime: " << delta_t << std::endl;
    }
}

TEST_CASE("Fluid_Surface", "") {
    SECTION("FluidSurface") {
        const double particleDiameter = 0.1;

        // Sample Particles in a Box
        FluidSystem particles =
            learnSPH::System::ParticleEmitter::getInstance().sampleFluidBox(
                Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(1, 1, 1),
                particleDiameter);
        const double r =
            learnSPH::Kernel::CubicSpline::support(particles.getSmoothingLength());
        std::shared_ptr<CompactNSearch::NeighborhoodSearch> nsearch;
        nsearch = std::make_shared<CompactNSearch::NeighborhoodSearch>(r);
        particles.addToNeighborhood(nsearch);
        nsearch->find_neighbors();
        particles.updateNormalizedDensities();
        std::vector<double> gridSDF;
        std::vector<Eigen::Vector3d> gridVerts;
        Eigen::Vector3i gridDims;
        discretizeFluidSystemSDF(
            particles.getPositions(), particles.getNormalizedDensities(),
            particles.getKernelLookUp(), particles.getSmoothingLength(), 0.6,
            particles.getSmoothingLength(), &gridSDF, &gridVerts, &gridDims);

        std::vector<Eigen::Vector3d> vertices;
        std::vector<std::array<int, 3>> triangles;
        marchCubes(gridDims, gridSDF, gridVerts, vertices, triangles);

        std::stringstream filename;
        filename << SOURCE_DIR << "/res/surface/fluid_surface.vtk";

        save_mesh_to_vtk(filename.str(), vertices, triangles);
    }
}

TEST_CASE("SimulationSurface", "") {
    SECTION("SurfaceSimulation") {
        double ratioSmoothingLengthSamplingStep = 2.0;
        const double particleDiameter = 0.05; // Reset to 0.05

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
                               Eigen::Vector3d(2.5, 2.5, 2.5), // Reset to 2.5, 2.5, 2.5
                               particleDiameter));
        solver.setBoundaryViscosity(0, 0.02);

        std::vector<SurfaceInformation> surfaceInfos;
        solver.run("complex_simulation", 6000, &surfaceInfos);
        for (SurfaceInformation surfaceInfo : surfaceInfos) {
            std::vector<double> gridSDF;
            std::vector<Eigen::Vector3d> gridVerts;
            Eigen::Vector3i gridDims;
            discretizeFluidSystemSDF(
                surfaceInfo.getPositions(), surfaceInfo.getNormalizedDensities(),
                surfaceInfo.getKernelLookup(), surfaceInfo.getSmoothingLength(), 0.6,
                surfaceInfo.getSmoothingLength() / ratioSmoothingLengthSamplingStep,
                &gridSDF, &gridVerts, &gridDims);

            std::vector<Eigen::Vector3d> vertices;
            std::vector<std::array<int, 3>> triangles;
            marchCubes(gridDims, gridSDF, gridVerts, vertices, triangles);
            std::stringstream absoluteFilename;
            absoluteFilename << SOURCE_DIR << "/res/simulation/"
                             << ratioSmoothingLengthSamplingStep
                             << surfaceInfo.getFilename() << ".vtk";
            save_mesh_to_vtk(absoluteFilename.str(), vertices, triangles);
        }
    }
}

