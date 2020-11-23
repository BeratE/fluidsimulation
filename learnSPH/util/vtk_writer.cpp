#include "vtk_writer.h"

void learnSPH::save_particles_to_vtk(std::string path,
                                     const std::vector<Eigen::Vector3d> &positions)
{
    save_particles_to_vtk(path, positions, std::vector<double>(), std::vector<Eigen::Vector3d>());
}


void learnSPH::save_particles_to_vtk(std::string path,
                           const std::vector<Eigen::Vector3d> &positions,
                           const std::vector<double> &scalar_data)
{
    save_particles_to_vtk(path, positions, scalar_data, std::vector<Eigen::Vector3d>());
}

void learnSPH::save_particles_to_vtk(std::string path,
                                     const std::vector<Eigen::Vector3d>& positions,
                                     const std::vector<double>& scalar_data,
                                     const std::vector<Eigen::Vector3d>& vector_data)
{
    if (positions.size() == 0)
        return;
        
    vtkio::VTKFile vtk_file;

        vtk_file.set_points_from_twice_indexable(positions);
        // Convenient way to specify that we don't have cells, just points
        vtk_file.set_cells_as_particles(positions.size());

        if (scalar_data.size() > 0)
            vtk_file.set_point_data_from_indexable("scalar",
                                                   scalar_data,
                                                   vtkio::AttributeType::Scalars);

        if (vector_data.size() > 0)
            vtk_file.set_point_data_from_twice_indexable("vector",
                                                         vector_data,
                                                         vtkio::AttributeType::Vectors);
    vtk_file.write(path);
}
