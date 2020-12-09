#if (defined __GNUC__) && (__GNUC__>4 || __GNUC_MINOR__>=7)
  #undef _GLIBCXX_ATOMIC_BUILTINS
  #undef _GLIBCXX_USE_INT128
#endif

#include "surface.h"
#include "mc_lut.h"
#include "helper.cu"
#include <iostream>


// Constant data
__constant__ dim3 VOL_DIM;


// Device Code
__inline__ __device__
uint getVertIdx(uint x, uint y, uint z)
{
    return x * VOL_DIM.y * VOL_DIM.z + y * VOL_DIM.z + z;
}

__inline__ __device__
uint getCellIdx(uint x, uint y, uint z)
{
    return x * (VOL_DIM.y -1) * (VOL_DIM.z -1) + y * (VOL_DIM.z -1) + z;
}

__inline __device__
uint3 getCellPos(uint idx, dim3 volDim)
{
    uint x = idx/ ((VOL_DIM.y -1) * (VOL_DIM.z -1));
    uint r = idx - (x * ((VOL_DIM.y -1) * (VOL_DIM.z -1)));
    uint y = r / (volDim.z -1);
    uint z = r - (y * (volDim.z -1));
    return make_uint3(x, y, z);
}


__global__
void triangulateCell(int* edges, int *lut,
                     float3 *volVerts, float *volSDF,
                     float3 *outTrianglePoints)
{
    uint cellIdx = getGlobalIdx1d1d();

    // Kill unused threads
    if (cellIdx > ((VOL_DIM.x -1) * (VOL_DIM.y -1) * (VOL_DIM.z -1)))
        return;
    
    dim3 pos = getCellPos(cellIdx, VOL_DIM);
    
    //printf("id: %d, x: %d, y: %d, z: %d\n", idx, pos.x, pos.y, pos.z);   
    
    // Get vertex indices of the cell;
    uint vertIdx[8];
    vertIdx[0] = getVertIdx(pos.x,   pos.y,   pos.z);
    vertIdx[1] = getVertIdx(pos.x+1, pos.y,   pos.z);
    vertIdx[2] = getVertIdx(pos.x+1, pos.y+1, pos.z);
    vertIdx[3] = getVertIdx(pos.x,   pos.y+1, pos.z);
    vertIdx[4] = getVertIdx(pos.x,   pos.y,   pos.z+1);
    vertIdx[5] = getVertIdx(pos.x+1, pos.y,   pos.z+1);
    vertIdx[6] = getVertIdx(pos.x+1, pos.y+1, pos.z+1);
    vertIdx[7] = getVertIdx(pos.x,   pos.y+1, pos.z+1);
    
    uint caseIdx = 0;
    for (uint i = 0; i < 8; i++) 
        caseIdx = caseIdx | (1 << i) * (volSDF[vertIdx[i]] <= 0.0);
    
    // Iterate edges                
    for (uint i = 0; i < 5; i++) { // triangle number
        for (uint j = 0; j < 3; j++) { // triangle point number
            uint offset = cellIdx*5*3 + i*3 + j;
            int edgeIdx = lut[caseIdx*5*3 + i*3 + j];
            
            if (edgeIdx == -1) {
                outTrianglePoints[offset]
                    = make_float3(offset, offset, offset);
                    continue;
            }

            int* edge = &edges[edgeIdx*2];
            
            //Global 1d vertex indices
            size_t vA = vertIdx[edge[0]];
            size_t vB = vertIdx[edge[1]];
            // Isolevel vales
            float isoA = volSDF[vA];
            float isoB = volSDF[vB];
            // Interpolate positions
            float alpha = isoA / (isoA - isoB);

            outTrianglePoints[offset] = lerp(volVerts[vA], volVerts[vB], alpha);
        }
    }
}

//
// Host code 
//

using namespace learnSPH;

std::vector<Eigen::Vector3d>
Surface::marchCubes(Eigen::Vector3i volDim,
                    std::vector<Eigen::Vector3f> &volVerts,
                    std::vector<float> &volSDF)
{
    assert(volVerts.size() == volSDF.size());

    const size_t VERT_COUNT = volVerts.size();
    const size_t CELL_COUNT = (volDim(0)-1)*(volDim(1)-1)*(volDim(2)-1);
    
    // Allocate device memory    
    float3 *d_volVerts;
    float *d_volSDF;
    int *d_lut;
    int *d_edges;
    cudaMalloc((void **)&d_volVerts, VERT_COUNT * sizeof(float3));
    cudaMalloc((void **)&d_volSDF, VERT_COUNT * sizeof(float));
    cudaMalloc((void **)&d_lut, 256*5*3 * sizeof(int));
    cudaMalloc((void **)&d_edges, 12*2 * sizeof(int));
    
    // Allocate unified memory
    float3* hd_outTriangles;
    cudaMallocManaged(&hd_outTriangles,
                      3 * 5 * CELL_COUNT * sizeof(float3));
    
    // Copy memory to device
    cudaMemcpy(d_volVerts, volVerts.front().data(),
               VERT_COUNT * sizeof(float3), cudaMemcpyHostToDevice);
    cudaMemcpy(d_volSDF, volSDF.data(),
               VERT_COUNT * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_lut, &MARCHING_CUBES_TABLE,
               256*5*3 * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_edges, &CELL_EDGES,
               12*2 * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(VOL_DIM, volDim.data(), sizeof(dim3));

    // Run kernel
    uint blockSize = 256;
    uint gridSize = (CELL_COUNT + (blockSize -1)) /blockSize;
    //uint sharedMemSize = VERT_COUNT * sizeof(float3);

#ifndef NDEBUG
    std::cout << "BlockSize: " << blockSize << "\nGridSize: " <<  gridSize << std::endl;
#endif
    
    triangulateCell<<<gridSize, blockSize>>>(
        d_edges, d_lut,
        d_volVerts, d_volSDF,
        hd_outTriangles);

    cudaDeviceSynchronize();
    
    // Reconstruct triangles
    std::vector<Eigen::Vector3d> triangles;
    for (uint c = 0; c < CELL_COUNT; c++) {
        for (uint i = 0; i < 5; i++) {
            for (uint j = 0; j < 3; j++) {
                uint offset = c*5*3 + i*3 + j;
                
                float3 point = hd_outTriangles[offset];
                                
                if (point.x == offset && point.y == offset && point.z == offset)
                    continue;

#ifndef NDEBUG
                std::cout << "N " << offset << ": "
                          << point.x << " " << point.y << " " << point.z << std::endl;
#endif                
                

                triangles.push_back(Eigen::Vector3d(point.x, point.y, point.z));
            }
        }
    }
    
#ifndef NDEBUG
    std::cout << "Number of Triangles: " << triangles.size() << std::endl;
#endif

    // Free device memory
    cudaFree(&d_volVerts);
    cudaFree(&d_volSDF);
    cudaFree(&d_lut);
    cudaFree(&d_edges);
    cudaFree(&hd_outTriangles);

    return triangles;
}
