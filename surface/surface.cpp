#include "surface.h"
#include <iostream>
#include <learnSPH/kernel.h>
#include <math.h>
#include <cfloat>

using namespace learnSPH;
using namespace learnSPH::Surface;

void Surface::discretizeSDF(
    const Eigen::Vector3d volSize,
    const Eigen::Vector3i volDims,
    std::function<double(Eigen::Vector3d)> const &sdf,
    std::vector<double> *pOutVolSDF,
    std::vector<Eigen::Vector3d> *pOutVolVerts)
{
    Eigen::Vector3d stepSizes = volSize;
    stepSizes(0) /= (double)volDims(0);
    stepSizes(1) /= (double)volDims(1);
    stepSizes(2) /= (double)volDims(2);

    for (size_t x = 0; x < volDims(0); x++) {
        for (size_t y = 0; y < volDims(1); y++) {
            for (size_t z = 0; z < volDims(2); z++) {
                pOutVolVerts->push_back(
                    Eigen::Vector3d(
                        stepSizes(0) * x,
                        stepSizes(1) * y,
                        stepSizes(0) * z));
                
                pOutVolSDF->push_back(sdf(pOutVolVerts->back()));
            }
        }
    }
}

void Surface::discretizeFluidSystemSDF(
    const std::vector<Eigen::Vector3d>& positions,
    const std::vector<double>& normalizedDensities,
    const Kernel::CubicSpline::Table& kernelLookup,
    double smoothingLength,
    double samplingDistance,
    double paramC,
    std::vector<double>* pOutVolSDF,
    std::vector<Eigen::Vector3d>* pOutVolVerts,
    Eigen::Vector3i* pOutDims)
{
    // TODO: Find Bounding Box for the fluid
    Eigen::Vector3d lowerLeft = Eigen::Vector3d::Zero();
    Eigen::Vector3d upperRight = Eigen::Vector3d::Zero();
    Surface::boundingBox(positions, lowerLeft, upperRight);

    // TODO: Increase the length of the bounding box by 3*samplingDistance in every direction
    lowerLeft -= 2 * smoothingLength * Eigen::Vector3d::Ones();
    upperRight += 2 * smoothingLength * Eigen::Vector3d::Ones();

    (*pOutDims) = Eigen::Vector3i(
        (upperRight.x() - lowerLeft.x()) / samplingDistance,
        (upperRight.y() - lowerLeft.y()) / samplingDistance,
        (upperRight.z() - lowerLeft.z()) / samplingDistance
        );

    // TODO: Create Vertices and initialize SDF values for bounding box
    for (size_t x = 0; x < pOutDims->x(); x++) {
        for (size_t y = 0; y < pOutDims->y(); y++) {
            for (size_t z = 0; z < pOutDims->z(); z++) {
                pOutVolSDF->push_back(-paramC);
                pOutVolVerts->push_back(
                    Eigen::Vector3d(
                        lowerLeft.x() + x * samplingDistance,
                        lowerLeft.y() + y * samplingDistance,
                        lowerLeft.z() + z * samplingDistance));
            }
        }
    }

    // TODO: Calculate number of steps required in every direction (positive and negative)
    double support = learnSPH::Kernel::CubicSpline::support(smoothingLength);
    int samplingSteps = ceil(support / samplingDistance);

    #pragma omp parallel for schedule(static)
    for (int posIdx = 0; posIdx < positions.size(); posIdx++) {
        const Eigen::Vector3d& position = positions[posIdx];
        // TODO: Find the indiizes that belongs to the lower left vertex of the cell that the particle is positioned in 
        Eigen::Vector3d offset = position - lowerLeft;
        Eigen::Vector3i indizes = Eigen::Vector3i(
            floor(offset.x() / samplingDistance),
            floor(offset.y() / samplingDistance),
            floor(offset.z() / samplingDistance));

        // TODO: Sample from the bounding box using the determined numbers of steps in each direction
        for (int xStepNr = -samplingSteps; xStepNr <= samplingSteps; xStepNr++) {
            for (int yStepNr = -samplingSteps; yStepNr <= samplingSteps; yStepNr++) {
                for (int zStepNr = -samplingSteps; zStepNr <= samplingSteps; zStepNr++) {
                    Eigen::Vector3i gridVertexIndizes = Eigen::Vector3i(
                        indizes.x() + xStepNr,
                        indizes.y() + yStepNr,
                        indizes.z() + zStepNr
                        );
                    if (gridVertexIndizes.x() < 0 || gridVertexIndizes.x() >= pOutDims->x() ||
                        gridVertexIndizes.y() < 0 || gridVertexIndizes.y() >= pOutDims->y() ||
                        gridVertexIndizes.z() < 0 || gridVertexIndizes.z() >= pOutDims->z()) {
                        continue;
                    }
                    size_t gridVertexIdx = Surface::getVertIdx(gridVertexIndizes, (*pOutDims));
                    Eigen::Vector3d gridVertex = pOutVolVerts->at(gridVertexIdx);
                    
                    if ((gridVertex - position).norm() < support) {
                        (*pOutVolSDF)[gridVertexIdx] += (1.0 / normalizedDensities[posIdx]) * kernelLookup.weight(gridVertex, position);
                    }
                }
            }
        }
    }
}

void Surface::boundingBox(const std::vector<Eigen::Vector3d>& positions, Eigen::Vector3d& bottomLeft, Eigen::Vector3d& upperRight) {
    double minX = DBL_MAX, minY = DBL_MAX, minZ = DBL_MAX;
    double maxX = DBL_MIN, maxY = DBL_MIN, maxZ = DBL_MIN;
    for (auto &position : positions) {
        if (position.x() < minX)
            minX = position.x();
        if (position.y() < minY)
            minY = position.y();
        if (position.z() < minZ)
            minZ = position.z();
        if (position.x() > maxX)
            maxX = position.x();
        if (position.y() > maxY)
            maxY = position.y();
        if (position.z() > maxZ)
            maxZ = position.z();
    }
    bottomLeft(0) = minX;
    bottomLeft(1) = minY;
    bottomLeft(2) = minZ;
    upperRight(0) = maxX;
    upperRight(1) = maxY;
    upperRight(2) = maxZ;
}
