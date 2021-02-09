#pragma once
#include <vector>
#include <array>
#include <fstream>
#include <iostream>

#include <Eigen/Dense>
#include <vtkio/VTKFile.h>

namespace learnSPH
{
    class Solver;
    const std::string lblDensities = "Densities";
    
    void readDensitiesFromVTK(
        std::string path,
        std::vector<double>& densities);

    void writeParticlesToVTK(
        std::string path,
        const std::vector<Eigen::Vector3d> &positions,
        const std::vector<double> &scalar_data,
        std::string comments = "");

    void writeMeshToVTK(
        std::string path,
        const std::vector<Eigen::Vector3d> &vertices,
        const std::vector<std::array<int, 3>> &triangles);

    void readParticlesFromVTK(
        std::string path,
        std::vector<Eigen::Vector3d> &positions,
        std::string &comments);

    void outputParams(std::string filename, Solver &solver, double runTime);
}
