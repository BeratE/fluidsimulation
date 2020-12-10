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


static void compute_node_normals(std::vector<Eigen::Vector3f>& out_normals,
                                  const std::vector<Eigen::Vector3f>& vertices,
                                  const std::vector<std::array<int, 3>>& triangles)
{
    out_normals.clear();
    out_normals.resize(vertices.size(), Eigen::Vector3f::Zero());
    for (const std::array<int, 3> &triangle : triangles) {
        const Eigen::Vector3f& p0 = vertices[triangle[0]];
        const Eigen::Vector3f& p1 = vertices[triangle[1]];
        const Eigen::Vector3f& p2 = vertices[triangle[2]];

        const Eigen::Vector3f triangle_normal = (p0 - p2).cross(p1 - p2);
        const float triangle_area =
            0.5 * std::abs((p1[0] - p0[0]) * (p2[1] - p0[1])
                           - (p2[0] - p0[0]) * (p1[1] - p0[1]));

        out_normals[triangle[0]] += triangle_area * triangle_normal;
        out_normals[triangle[1]] += triangle_area * triangle_normal;
        out_normals[triangle[2]] += triangle_area * triangle_normal;
    }

    for (Eigen::Vector3f& normal : out_normals) {
        normal.normalize();
    }
}



void learnSPH::save_mesh_to_vtk(std::string path,
                                const std::vector<Eigen::Vector3f>& vertices,
                                const std::vector<std::array<int, 3>>& triangles)
{
    vtkio::VTKFile vtk_file;

    vtk_file.set_points_from_twice_indexable(vertices);
    vtk_file.set_cells_from_twice_indexable(triangles, vtkio::CellType::Triangle);

    std::vector<Eigen::Vector3f> normals;
    compute_node_normals(normals, vertices, triangles);
    vtk_file.set_point_data_from_twice_indexable("normals", normals, vtkio::AttributeType::Normals);

    vtk_file.write(path);
}
