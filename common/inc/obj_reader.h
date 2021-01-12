#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <array>


#include <Eigen/Dense>

namespace learnSPH
{
    struct TriMesh
    {
        std::vector<Eigen::Vector3d> vertices;
        std::vector<std::array<int, 3>> triangles;
    };

    std::vector<TriMesh> readTriMeshesFromObj(std::string filename);
}
