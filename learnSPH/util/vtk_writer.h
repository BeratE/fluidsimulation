#pragma once
#include <vector>
#include <array>
#include <fstream>
#include <iostream>

#include <Eigen/Dense>
#include <vtkio/VTKFile.h>

namespace learnSPH
{
	void save_particles_to_vtk(std::string path, const std::vector<Eigen::Vector3d>& particles, const std::vector<double>& particle_scalar_data, const std::vector<Eigen::Vector3d>& particle_vector_data);
}
