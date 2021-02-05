#pragma once
#include <vector>
#include <array>
#include <fstream>
#include <iostream>

#include <Eigen/Dense>
#include <vtkio/VTKFile.h>

namespace learnSPH
{
    const std::string lblDensities = "Densities";
    
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
}
