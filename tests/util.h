#pragma once
#include "learnSPH/solver.h"
#include "learnSPH/system/particlesystem.h"
#include "learnSPH/system/boundarysystem.h"
#include <iostream>

using namespace learnSPH;
using namespace learnSPH::System;

void outputParams(std::string filename, Solver &solver, double runTime)
{
    std::stringstream output;
    output << filename << std::endl;
    output << "Num Particles: " << solver.getSystem().getSize() << std::endl;
    output << "Runtime: " << runTime << std::endl;
    output << "Max Timestep: " << solver.getMaxTimeStepSeconds() << std::endl;
    output << "XSPH Smoothing: " << solver.getParamSmoothing() << std::endl;
    output << "Smoothing Enabled: " << solver.smoothingEnabled() << std::endl;
    output << "Gravity Enabled: " << solver.gravityEnabled() << std::endl;
    output << "Tension Enabled: " << solver.tensionEnabled() << std::endl;
    output << "Adhesion Enabled: " << solver.adhesionEnabled() << std::endl;
    output << "Drag Enabled: " << solver.dragEnabled() << std::endl;
    output << "Drag : " << solver.getParamDrag() << std::endl;

    output << "Fluid Rest Density: " << solver.getSystem().getRestDensity() << std::endl;
    output << "Fluid Surface Tension: " << solver.getSystem().getGamma() << std::endl;
    output << "Fluid Viscosity: " << solver.getSystem().getViscosity() << std::endl;

    size_t i = 0;
    for (auto const &boundary : solver.getBoundaries()) {
        output << "Boundary " << i << " Rest Density: "  << boundary.getRestDensity() << std::endl;
        output << "Boundary " << i << " Viscosity: "  << boundary.getViscosity() << std::endl;
        output << "Boundary " << i << " Adhesion: " << boundary.getBeta() << std::endl;
        i++;
    }

    std::stringstream name;
    name << filename << ".out";
    
    std::ofstream outfile;
    outfile.open(name.str());
    outfile << output.rdbuf();
    outfile.close();

    std::cout << output.str() << std::endl;
}
