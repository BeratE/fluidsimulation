#pragma once
#include <vector>
#include <array>
#include <fstream>
#include <iostream>

#include <Eigen/Dense>
#include <vtkio/VTKFile.h>

namespace learnSPH
{
    void save_particles_to_vtk(std::string path,
                               const std::vector<Eigen::Vector3d>& positions);
    void save_particles_to_vtk(std::string path,
                               const std::vector<Eigen::Vector3d>& positions,
                               const std::vector<double>& scalar_data);
    void save_particles_to_vtk(std::string path,
                               const std::vector<Eigen::Vector3d>& positions,
                               const std::vector<double>& scalar_data,
                               const std::vector<Eigen::Vector3d>& vector_data);

    void save_mesh_to_vtk(std::string path,
                          const std::vector<Eigen::Vector3d>& vertices,
                          const std::vector<std::array<int, 3>>& triangles);
}
