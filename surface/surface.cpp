#include "surface.h"
#include <iostream>
#include <learnSPH/kernel.h>
#include <math.h>

using namespace learnSPH;

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
    const learnSPH::System::FluidSystem& system,
    const double c,
    const double samplingDistance,
    std::vector<double>* pOutVolSDF,
    std::vector<Eigen::Vector3d>* pOutVolVerts,
    Eigen::Vector3i* pOutDims) {
    // TODO: Find Bounding Box for the fluid
    Eigen::Vector3d lowerLeft = Eigen::Vector3d::Zero();
    Eigen::Vector3d upperRight = Eigen::Vector3d::Zero();
    system.boundingBox(lowerLeft, upperRight);

    // TODO: Increase the length of the bounding box by 3*samplingDistance in every direction
    lowerLeft -= 2 * samplingDistance * Eigen::Vector3d::Ones();
    upperRight += 2 * samplingDistance * Eigen::Vector3d::Ones();

    (*pOutDims) = Eigen::Vector3i(
        (upperRight.x() - lowerLeft.x()) / samplingDistance,
        (upperRight.y() - lowerLeft.y()) / samplingDistance,
        (upperRight.z() - lowerLeft.z()) / samplingDistance
    );

    // TODO: Create Vertices and initialize SDF values for bounding box
    for (size_t x = 0; x < pOutDims->x(); x++) {
        for (size_t y = 0; y < pOutDims->y(); y++) {
            for (size_t z = 0; z < pOutDims->z(); z++) {
                pOutVolSDF->push_back(-c);
                pOutVolVerts->push_back(Eigen::Vector3d(
                    lowerLeft.x() + x * samplingDistance,
                    lowerLeft.y() + y * samplingDistance,
                    lowerLeft.z() + z * samplingDistance
                ));
            }
        }
    }


    // TODO: Calculate number of steps required in every direction (positive and negative)
    double support = learnSPH::Kernel::CubicSpline::support(system.getSmoothingLength());
    int samplingSteps = ceil(support / samplingDistance);
    std::cout << samplingSteps << std::endl;
    for (size_t posIdx = 0; posIdx < system.getSize(); posIdx++) {
        const Eigen::Vector3d& position = system.getPositions()[posIdx];
        // TODO: Find the indiizes that belongs to the lower left vertex of the cell that the particle is positioned in 
        Eigen::Vector3d offset = position - lowerLeft;
        Eigen::Vector3i indizes = Eigen::Vector3i(
            floor(offset.x() / samplingDistance),
            floor(offset.y() / samplingDistance),
            floor(offset.z() / samplingDistance)
        );

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
                    size_t gridVertexIdx = getVertIdx(gridVertexIndizes, (*pOutDims));
                    Eigen::Vector3d gridVertex = pOutVolVerts->at(gridVertexIdx);
                    
                    if ((gridVertex - position).norm() < support) {
                        (*pOutVolSDF)[gridVertexIdx] += (1.0 / system.getNormalizedDensities()[posIdx]) * system.calculateWeightBetweenParticles(gridVertex, position);
                    }
                }
            }
        }
    }
}

