#include "surface.h"
#include <iostream>

using namespace learnSPH;

void Surface::discretizeSDF(
    const Eigen::Vector3f volSize,
    const Eigen::Vector3i volDims,
    std::function<float(Eigen::Vector3f)> const &sdf,
    std::vector<float> *pOutVolSDF,
    std::vector<Eigen::Vector3f> *pOutVolVerts)
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
