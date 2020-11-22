#include "vtk_writer.h"

void learnSPH::save_particles_to_vtk(std::string path,
                           const std::vector<Eigen::Vector3d> &particles,
                           const std::vector<double> &particle_scalar_data)
{
    save_particles_to_vtk(path, particles, particle_scalar_data, std::vector<Eigen::Vector3d>());
}

void learnSPH::save_particles_to_vtk(std::string path,
                                     const std::vector<Eigen::Vector3d>& particles,
                                     const std::vector<double>& particle_scalar_data,
                                     const std::vector<Eigen::Vector3d>& particle_vector_data)
{
    vtkio::VTKFile vtk_file;
    vtk_file.set_points_from_twice_indexable(particles);
    // Convenient way to specify that we don't have cells, just points
    vtk_file.set_cells_as_particles(particles.size());
    
    if (particle_scalar_data.size() > 0)
        vtk_file.set_point_data_from_indexable("scalar", particle_scalar_data, vtkio::AttributeType::Scalars);
    
    if (particle_vector_data.size() > 0)
        vtk_file.set_point_data_from_twice_indexable("vector", particle_vector_data, vtkio::AttributeType::Vectors);
    
    vtk_file.write(path);
}
