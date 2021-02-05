#include "vtk_writer.h"

// void VTK::writeParticlesToVTK(
//     std::string path,
//     const std::vector<Eigen::Vector3d> &positions)
// {
//     VTK::writeParticlesToVTK(path, positions, std::vector<double>(), std::vector<Eigen::Vector3d>());
// }

void learnSPH::writeParticlesToVTK(
    std::string path,
    const std::vector<Eigen::Vector3d> &positions,
    const std::vector<double> &scalar_data,
    std::string comments)
{
    if (positions.size() == 0)
        return;
        
    vtkio::VTKFile vtk_file;

    vtk_file.set_points_from_twice_indexable(positions);
    
    // Convenient way to specify that we don't have cells, just points
    vtk_file.set_cells_as_particles(positions.size());

    vtk_file.set_point_data_from_indexable("scalar", scalar_data,
                                           vtkio::AttributeType::Scalars);
    
    vtk_file.set_comments(comments);
    
    vtk_file.write(path);

    std::cout << "Saved results to " << path << std::endl;
}


static void compute_node_normals(
    std::vector<Eigen::Vector3d>& out_normals,
    const std::vector<Eigen::Vector3d>& vertices,
    const std::vector<std::array<int, 3>>& triangles)
{
    out_normals.clear();
    out_normals.resize(vertices.size(), Eigen::Vector3d::Zero());
    for (const std::array<int, 3> &triangle : triangles) {
        const Eigen::Vector3d& p0 = vertices[triangle[0]];
        const Eigen::Vector3d& p1 = vertices[triangle[1]];
        const Eigen::Vector3d& p2 = vertices[triangle[2]];

        const Eigen::Vector3d triangle_normal = (p0 - p2).cross(p1 - p2);
        const double triangle_area =
            0.5 * std::abs((p1[0] - p0[0]) * (p2[1] - p0[1])
                           - (p2[0] - p0[0]) * (p1[1] - p0[1]));

        out_normals[triangle[0]] += triangle_area * triangle_normal;
        out_normals[triangle[1]] += triangle_area * triangle_normal;
        out_normals[triangle[2]] += triangle_area * triangle_normal;
    }

    for (Eigen::Vector3d& normal : out_normals) {
        normal.normalize();
    }
}

void learnSPH::writeMeshToVTK (
    std::string path,
    const std::vector<Eigen::Vector3d>& vertices,
    const std::vector<std::array<int, 3>>& triangles)
{
    vtkio::VTKFile vtk_file;

    vtk_file.set_points_from_twice_indexable(vertices);
    vtk_file.set_cells_from_twice_indexable(triangles, vtkio::CellType::Triangle);

    std::vector<Eigen::Vector3d> normals;
    compute_node_normals(normals, vertices, triangles);
    vtk_file.set_point_data_from_twice_indexable("normals", normals, vtkio::AttributeType::Normals);

    vtk_file.write(path);

    std::cout << "Saved results to " << path << std::endl;
}

void learnSPH::readParticlesFromVTK(
    std::string path,
    std::vector<Eigen::Vector3d> &positions,
    std::string &comments)
{
    vtkio::VTKFile vtk_file;
    
    vtk_file.read(path);
   
    const size_t numPoints = vtk_file.get_number_of_points();
    positions.resize(numPoints);    

    vtk_file.get_points_to_twice_indexable(positions);
    comments = vtk_file.get_comments();
}
