#include "surface.h"
#include <iostream>

using namespace learnSPH;

void Surface::discretizeSDF(
    const Eigen::Vector3f volSize,
    const Eigen::Vector3i volDims,
    std::function<float(Eigen::Vector3f)> const &sdf,
    std::vector<Eigen::Vector3f> *pOutVolVerts,
    std::vector<float> *pOutVolSDF)
{
    Eigen::Vector3f stepSizes = volSize;
    stepSizes(0) /= (float)volDims(0);
    stepSizes(1) /= (float)volDims(1);
    stepSizes(2) /= (float)volDims(2);

    for (size_t x = 0; x < volDims(0); x++) {
        for (size_t y = 0; y < volDims(1); y++) {
            for (size_t z = 0; z < volDims(2); z++) {
                pOutVolVerts->push_back(
                    Eigen::Vector3f(
                        stepSizes(0) * x,
                        stepSizes(1) * y,
                        stepSizes(0) * z));
                
                pOutVolSDF->push_back(sdf(pOutVolVerts->back()));
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
