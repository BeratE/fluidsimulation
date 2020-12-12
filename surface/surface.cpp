#include "surface.h"
#include <iostream>
#include <learnSPH/kernel.h>
#include <math.h>

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

void Surface::discretizeFluidSystemSDF(
    const learnSPH::System::FluidSystem& system,
    const double c,
    const double samplingDistance,
    std::vector<double>* pOutVolSDF,
    std::vector<Eigen::Vector3d>* pOutVolVerts,
    CompactNSearch::NeighborhoodSearch& nsearch) {

    // TODO: Find Bounding Box for the fluid
    Eigen::Vector3d lowerLeft = Eigen::Vector3d::Zero();
    Eigen::Vector3d upperRight = Eigen::Vector3d::Zero();
    system.boundingBox(lowerLeft, upperRight);

    // TODO: Increase the length of the bounding box by a bit 3*samplingDistance in every direction
    lowerLeft -= 1.5 * samplingDistance * Eigen::Vector3d::Ones();
    upperRight += 1.5 * samplingDistance * Eigen::Vector3d::Ones();

    Eigen::Vector3i stepNrMax(
        (upperRight.x() - lowerLeft.x()) / samplingDistance,
        (upperRight.y() - lowerLeft.y()) / samplingDistance,
        (upperRight.z() - lowerLeft.z()) / samplingDistance
    );


    // TODO: Create Vertices and initialize SDF values for bounding box
    for (size_t x = 0; x < stepNrMax.x(); x++) {
        for (size_t y = 0; y < stepNrMax.y(); y++) {
            for (size_t z = 0; z < stepNrMax.z(); z++) {
                pOutVolSDF->push_back(-c);
                pOutVolVerts->push_back(Eigen::Vector3d(
                    lowerLeft.x() + x * samplingDistance,
                    lowerLeft.y() + y * samplingDistance,
                    lowerLeft.z() + z * samplingDistance
                ));
            }
        }
    }

    // TODO: Add Grid to NeighborHoodSearch
    unsigned int gridPSid = nsearch.add_point_set(pOutVolVerts->front().data(), pOutVolVerts->front().size());
    nsearch.find_neighbors();

    std::vector<double> normalizedDensities;
    calculateNormalizedDensities(system, nsearch, normalizedDensities);

    // TODO: Calculate number of steps required in every direction (positive and negative)
    double support = learnSPH::Kernel::CubicSpline::support(system.getSmoothingLength());
    int samplingSteps = floor(support / samplingDistance);

    for (const Eigen::Vector3d& position : system.getPositions()) {
        // TODO: Find the indiizes that belongs to the lower left vertex of the cell that the particle is positioned in 
        Eigen::Vector3d offset = position - lowerLeft;
        Eigen::Vector3i indizes = Eigen::Vector3i(
            floor(offset.x() / samplingDistance),
            floor(offset.y() / samplingDistance),
            floor(offset.z() / samplingDistance)
        );
        // TODO: Sample from the bounding box using the determined numbers of steps in each direction
        for (size_t xStepNr = -samplingSteps; xStepNr <= samplingSteps; xStepNr++) {
            for (size_t yStepNr = -samplingSteps; yStepNr <= samplingSteps; yStepNr++) {
                for (size_t zStepNr = -samplingSteps; zStepNr <= samplingSteps; zStepNr++) {
                    Eigen::Vector3i gridVertexIndizes = Eigen::Vector3i(
                        indizes.x() + xStepNr,
                        indizes.y() + yStepNr,
                        indizes.z() + zStepNr
                    );
                    if (gridVertexIndizes.x() < 0 || gridVertexIndizes.x() >= stepNrMax.x() ||
                        gridVertexIndizes.y() < 0 || gridVertexIndizes.y() >= stepNrMax.y() ||
                        gridVertexIndizes.z() < 0 || gridVertexIndizes.z() >= stepNrMax.z()) {
                        continue;
                    }
                    size_t gridVertexIdx = getVertIdx(gridVertexIndizes, stepNrMax);
                    Eigen::Vector3d gridVertex = pOutVolVerts->at(gridVertexIdx);
                    
                    if ((gridVertex - position).norm() < support) {
                        CompactNSearch::PointSet const& gridPS = nsearch.point_set(gridPSid);
                        for (size_t nIdx = 0; nIdx < gridPS.n_neighbors(system.getPointSetID(), gridVertexIdx); nIdx++) {
                            const unsigned int nId = gridPS.neighbor(system.getPointSetID(), gridVertexIdx, nIdx);
                            (*pOutVolSDF)[gridVertexIdx] += (1.0 / normalizedDensities[nId]) * system.calculateWeightBetweenParticles(gridVertex, system.getPositions()[nId]);
                        }
                    }
                }
            }
        }
    }
}

void Surface::calculateNormalizedDensities(const learnSPH::System::FluidSystem& system, CompactNSearch::NeighborhoodSearch& nsearch, std::vector<double>& normalizedDensities) {
    CompactNSearch::PointSet const ps = nsearch.point_set(system.getPointSetID());
    for (int fpId = 0; fpId < ps.n_points(); ++fpId) {
        double normalizedDensity = 0.0;
        for (size_t fnId = 0; ps.n_neighbors(system.getPointSetID(), fpId); fnId++) {
            const unsigned int nId = ps.neighbor(system.getPointSetID(), fpId, fnId);
            normalizedDensity += system.calculateWeightBetweenParticles(system.getPositions()[fpId], system.getPositions()[nId]);
        }
        normalizedDensities.push_back(normalizedDensity);
    }
}
