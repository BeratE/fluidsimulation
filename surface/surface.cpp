#include "surface.h"
#include "mc_lut.hpp"
#include <iostream>
#include <unordered_map>

using namespace learnSPH;

// void Surface::printCudaVersion()
// {
//     std::cout << "CUDA Compiled version: "
//               << __CUDACC_VER_MAJOR__ << "."
//               << __CUDACC_VER_MINOR__
//               << std::endl;

//     int runtime_ver;
//     cudaRuntimeGetVersion(&runtime_ver);
//     std::cout << "CUDA Runtime version: " << runtime_ver << std::endl;

//     int driver_ver;
//     cudaDriverGetVersion(&driver_ver);
//     std::cout << "CUDA Driver version: " << driver_ver << std::endl;
// }


inline size_t
getIndex1D(IndexVector3d it, Dimensions3d gridDim)
{
    return it(0) * gridDim(1) * gridDim(2) + it(1) * gridDim(2) + it(2);
}

std::vector<Eigen::Vector3d>
Surface::marchCubes(std::vector<Eigen::Vector3d> const &gridVerts,
                    std::vector<double> const &gridSDF,
                    const Dimensions3d gridDim)
{
    assert(gridSDF.size() == gridVerts.size());

    std::vector<Eigen::Vector3d> triangles;

    // Iterate cubes
    // Cube dimensions = number of vertices - 1
    for (size_t x = 0; x < gridDim(0)-1; x++) {
        for (size_t y = 0; y < gridDim(1)-1; y++) {
            for (size_t z = 0; z < gridDim(2)-1; z++) {
                // Mapping local index to global 1D grid index
                const std::array<size_t, 8> GRID_INDICES = {
                    getIndex1D(IndexVector3d(x, y, z), gridDim),
                    getIndex1D(IndexVector3d(x + 1, y, z), gridDim),
                    getIndex1D(IndexVector3d(x + 1, y + 1, z), gridDim),
                    getIndex1D(IndexVector3d(x, y + 1, z), gridDim),
                    getIndex1D(IndexVector3d(x, y, z + 1), gridDim),
                    getIndex1D(IndexVector3d(x + 1, y, z + 1), gridDim),
                    getIndex1D(IndexVector3d(x + 1, y + 1, z + 1), gridDim),
                    getIndex1D(IndexVector3d(x, y + 1, z + 1), gridDim),
                };
                               
                std::array<bool, 8> signs;
                for (size_t i = 0; i < 8; i++)
                    signs[i] = (gridSDF[GRID_INDICES[i]] <= 0.0);

                const std::array<std::array<int, 3>, 5> TRIANG
                    = getMarchingCubesCellTriangulation(signs);

                // Iterate edges                
                for (size_t i = 0; i < 5; i++) {
                    for (size_t j = 0; j < 3; j++) {
                        const int EDGE_INDEX = TRIANG[i][j];
                        if (EDGE_INDEX == -1)
                            break;                      

                        const auto EDGE = CELL_EDGES[EDGE_INDEX];
                        // Global 1d vertex indices
                        const size_t V_A = GRID_INDICES[EDGE[0]];
                        const size_t V_B = GRID_INDICES[EDGE[1]];
                        // Isolevel vales
                        const double ISO_A = gridSDF[V_A];
                        const double ISO_B = gridSDF[V_B];
                        // Interpolate positions
                        const double ALPHA = ISO_A / (ISO_A - ISO_B);
                        Eigen::Vector3d xs =
                            (1.0 - ALPHA) * gridVerts[V_A] + ALPHA * gridVerts[V_B];

                        triangles.push_back(xs);
                    }
                }                                
            }            
        }
    }

    return triangles;
}

void Surface::discretizeSDF(
    const Dimensions3d gridDim,
    const Eigen::Vector3d gridSize,
    std::function<double(Eigen::Vector3d)> const &sdf,
    std::vector<Eigen::Vector3d> *pOutGridVerts,
    std::vector<double> *pOutGridSDF)
{
    Eigen::Vector3d stepSizes = gridSize;
    stepSizes(0) /= (double)gridDim(0);
    stepSizes(1) /= (double)gridDim(1);
    stepSizes(2) /= (double)gridDim(2);

    for (size_t x = 0; x < gridDim(0); x++) {
        std::vector<Eigen::Vector3d> gridVerts;
        std::vector<double> gridSDF;
        for (size_t y = 0; y < gridDim(1); y++) {
            for (size_t z = 0; z < gridDim(2); z++) {
                pOutGridVerts->push_back(
                    Eigen::Vector3d(
                        stepSizes(0) * x,
                        stepSizes(1) * y,
                        stepSizes(0) * z));
                
                pOutGridSDF->push_back(sdf(pOutGridVerts->back()));
            }
        }
    }
}

// void constructSurface(const std::vector<Eigen::Vector3d> &gridVerts,
//                       const std::vector<double> &gridSDF,
//                       const Dimensions3d gridDim)
// {
//     assert(gridSDF.size() == gridVertices.size());

//     // Map Edge intex to intersection vertex
//     std::unordered_map<size_t, size_t> edgeIntersecs;
//     std::vector<Eigen::Vector3d> intersecVerts;

//     // Iterate Grid and look for intersections
//     for (size_t x = 0; x < gridDim(0); x++) {
//         for (size_t y = 0; y < gridDim(1); y++) {
//             for (size_t z = 0; z < gridDim(2); z++) {
//                 // Vertex origin
//                 IndexVector3d it(x, y, z); // index
//                 const size_t Va = getIndex1D(it, gridDim);
//                 const double sdfVa = gridSDF[Va];

//                 size_t E = (3 * Va) - 1; // Edge index
                
//                 // Loop Edges 0, 1, 2 (x, y, z) 
//                 for (i = 0; i < 3; i++) {
//                     // skip index if component is on the boundary
//                     if (it(i) >= (gridDim(i) - 1))
//                         continue;
                    
//                     E += 1;
//                     it(i) += 1;

//                     const size_t Vb = getIndex(it, gridDim);
//                     const double sdfVb = gridSDF[Va];

//                     // different signs
//                     if ((sdfVa / sdfVb) < 0.0) { 
//                         const double alpha = sdfVa / (sdfVa-sdfVb);
//                         const Eigen::Vector3d xs =
//                             (1.0 - alpha)*gridVerts[Va]
//                             + alpha      *gridVerts[Vb];
//                         intersecVerts.push_back(xs);
//                         edgeIntersecs.insert(
//                             std::make_pair<size_t, size_t>(
//                                 E, intersecVerts.size()-1));
//                     }
                    
//                     it(i) -= 1;
//                 }
//             }
//         }
//     }
// }
