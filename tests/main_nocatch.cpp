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

int main(int argc, char* argv[]) {
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
    filename << SOURCE_DIR << "/res/smallDam/" << "smallDam";

    filename << "_PBF";
    solver = new SolverPBF(particles);
    solver->addBoundary(box);

    solver->enableGravity(true);
    solver->enableAdhesion(true);
    solver->enableDrag(true);
    solver->enableSmoothing(true);
    solver->enableTension(true);

    solver->setParamDrag(0.025);
    filename << "_I";
    ((SolverPBF*)solver)->setNumIterations(7);
    solver->setMaxTimeStepSeconds(0.004);
    solver->setParamSmoothing(0.1);

    solver->setFluidViscosity(0.001);
    solver->enableTension(false);

    solver->setBoundaryViscosity(0, 0.001); // box
    solver->enableAdhesion(false);

    solver->setSnapShotAfterMS(1000.0 / 40);

    const double startTime = omp_get_wtime();

    solver->run(filename.str(), 4000);

    double endTime = omp_get_wtime();

    std::cout << "Runtime: " << endTime - startTime << std::endl;

    outputParams(filename.str(), *solver, endTime - startTime);

    delete solver;
}