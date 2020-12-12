#include "surface.h"
#include "mc_lut.h"
#include <unordered_map>
#include <iostream>
#include <array>

using namespace learnSPH;


size_t Surface::getVertIdx(Eigen::Vector3i pos, Eigen::Vector3i volDim)
{
    return pos(0) * volDim(1) * volDim(2)
        + pos(1) * volDim(2)
        + pos(2);
}

static inline
Eigen::Vector3i getPos(size_t idx, Eigen::Vector3i volDim)
{
    uint x = idx/ (volDim(1) * volDim(2));
    uint r = idx - (x * (volDim(1) * volDim(2)));
    uint y = r / volDim(2);
    uint z = r - (y * volDim(2));
    return Eigen::Vector3i(x, y, z);
}

void collectVertices(const Eigen::Vector3i volDim,
                     const std::vector<float> &volSDF,
                     const std::vector<Eigen::Vector3f> &volVerts,
                     std::vector<Eigen::Vector3f> &outVertices,
                     std::unordered_map<size_t, size_t> &outEdgeIdxToVertIdx)
{
    const size_t SIZE = volDim(0)*volDim(1)*volDim(2);
    for (size_t idx = 0; idx < SIZE; idx++) {
        Eigen::Vector3i originPos = getPos(idx, volDim);
        size_t originIdx = idx;
        float  originLvl = volSDF[idx];
                
        // Loop Edges 0, 1, 2 (x, y, z)
        for (size_t i = 0; i < 3; i++) {
            if (originPos(i) == volDim(i)-1)
                continue;
            
            Eigen::Vector3i oppositePos = originPos;
            oppositePos(i) += 1;
            size_t oppositeIdx = Surface::getVertIdx(oppositePos, volDim);
            float oppositeLvl = volSDF[oppositeIdx];


            if (originLvl != 0 && (oppositeLvl / originLvl) < 0.0) {
                size_t edgeIdx = 3*originIdx + i;
                float alpha = originLvl / (originLvl -oppositeLvl);
                Eigen::Vector3f vert =
                    (1.0 -alpha) * volVerts[originIdx]
                    + alpha * volVerts[oppositeIdx];
                
                outVertices.push_back(vert);
                outEdgeIdxToVertIdx.insert(std::make_pair(edgeIdx, outVertices.size()-1));
            }
        }
    }
}

void Surface::marchCubes(const Eigen::Vector3i volDim,
                         const std::vector<float> &volSDF,
                         const std::vector<Eigen::Vector3f> &volVerts,
                         std::vector<Eigen::Vector3f> &outVertices,
                         std::vector<std::array<int, 3>> &outTriangles)
{
    size_t SIZE = volDim(0)*volDim(1)*volDim(2);
    assert(SIZE == volSDF.size());
    assert(volSDF.size() == volVerts.size());    
   
    std::unordered_map<size_t, size_t> edgeIdxToVertIdx;    
    collectVertices(volDim, volSDF, volVerts, outVertices, edgeIdxToVertIdx);

    // for (auto &v : edgeIdxToVertIdx)
    //     std::cout << v.first << " " << v.second << std::endl;

    // Iterate cubes
    // Cube dimensions = number of vertices - 1
    Eigen::Vector3i numCells = volDim - Eigen::Vector3i(1, 1, 1);
    const size_t NUM_CELLS = numCells(0) * numCells(1) * numCells(2);
    for (size_t idx = 0; idx < NUM_CELLS; idx++) {
        Eigen::Vector3i pos = getPos(idx, numCells);

        // Mapping local index to global 1D grid index
        std::array<size_t, 8> vertIdx;
        for (int i = 0; i < 8; i++) {
            const Eigen::Vector3i offset(CELL_VERTICES[i]);
            vertIdx[i] = getVertIdx(pos + offset, volDim);
        }

        uint caseIdx = 0;
        for (uint i = 0; i < 8; i++) {
            if (volSDF[vertIdx[i]] <= 0.0)
                caseIdx |= (1 << i);
        }

        // Iterate edges        
        for (uint i = 0; i < 5; i++) {
            std::array<int, 3> triangle;
            for (uint j = 0; j < 3; j++) {
                uint edgeIdx = MARCHING_CUBES_TABLE[caseIdx][i][j];
                if (edgeIdx == -1)
                    goto endloop;
                
                size_t originIdx  = vertIdx[CELL_EDGES[edgeIdx][0]];
                size_t globalEdgeIdx = 3*originIdx + EDGE_DIR[edgeIdx];
                printf("%d", globalEdgeIdx);
                triangle[j] = edgeIdxToVertIdx[globalEdgeIdx];
                printf(", %d\n", triangle[j]);
            }
            outTriangles.push_back(triangle);
        }

      endloop:
        continue;
    }

    // for (auto &v : outTriangles) {
    //     printf("%d, %d, %d\n", v[0], v[1], v[2]);
    // }
}
