#include "surface.h"
#include "mc_lut.h"

using namespace learnSPH;

static inline
size_t getVertIdx(Eigen::Vector3i pos, Eigen::Vector3i volDims)
{
    return pos(0) * volDims(1) * volDims(2)
        + pos(1) * volDims(2)
        + pos(2);
}

std::vector<Eigen::Vector3d>
Surface::marchCubes(Eigen::Vector3i volDim,
                    std::vector<Eigen::Vector3f> &volVerts,
                    std::vector<float> &volSDF)
{
    assert(volSDF.size() == volVerts.size());

    std::vector<Eigen::Vector3d> triangles;

    // Iterate cubes
    // Cube dimensions = number of vertices - 1
    for (size_t x = 0; x < volDim(0)-1; x++) {
        for (size_t y = 0; y < volDim(1)-1; y++) {
            for (size_t z = 0; z < volDim(2)-1; z++) {
                // Mapping local index to global 1D grid index
                const std::array<size_t, 8> VERT_IDX = {
                    getVertIdx(Eigen::Vector3i(x,   y,   z), volDim),
                    getVertIdx(Eigen::Vector3i(x+1, y,   z), volDim),
                    getVertIdx(Eigen::Vector3i(x+1, y+1, z), volDim),
                    getVertIdx(Eigen::Vector3i(x,   y+1, z), volDim),
                    getVertIdx(Eigen::Vector3i(x,   y,   z+1), volDim),
                    getVertIdx(Eigen::Vector3i(x+1, y,   z+1), volDim),
                    getVertIdx(Eigen::Vector3i(x+1, y+1, z+1), volDim),
                    getVertIdx(Eigen::Vector3i(x,   y+1, z+1), volDim),
                };
                               
                uint index = 0;
                for (uint i = 0; i < 8; i++) {
                    if (volSDF[VERT_IDX[i]] <= 0.0)
                        index |= (1 << i);
                }

                // Iterate edges                
                for (uint i = 0; i < 5; i++) {
                    for (uint j = 0; j < 3; j++) {
                        uint edgeIdx = MARCHING_CUBES_TABLE[index][i][j];
                        if (edgeIdx == -1)
                            goto endloop;                      

                        // Global 1d vertex indices
                        const size_t V_A = VERT_IDX[CELL_EDGES[edgeIdx][0]];
                        const size_t V_B = VERT_IDX[CELL_EDGES[edgeIdx][1]];
                        // Isolevel vales
                        const float ISO_A = volSDF[V_A];
                        const float ISO_B = volSDF[V_B];
                        // Interpolate positions
                        const float ALPHA = ISO_A / (ISO_A - ISO_B);
                        Eigen::Vector3d xs =
                            (1.0 - ALPHA) * volVerts[V_A].cast<double>()
                            + ALPHA * volVerts[V_B].cast<double>();

                        triangles.push_back(xs);
                    }
                }

            endloop:
                continue;
            }
        }
    }

    return triangles;
}
