#include "emitter.h"

using namespace learnSPH;

FluidSystem Emitter::sampleFluidBox(Eigen::Vector3d bottomLeft,
                                    Eigen::Vector3d topRight,
                                    double fluidRestDensity,
                                    double samplingDistance)
{
    const Eigen::Vector3d diagonal = topRight - bottomLeft;
    const Eigen::Vector3d diagonalSign = diagonal.array() / diagonal.array().abs();
    const Eigen::Vector3d samplingDir = samplingDistance * diagonalSign;
    const Eigen::Vector3d numSamples = Eigen::floor(Eigen::abs(diagonal.array())
                                                    / samplingDistance);
    
    FluidSystem particles;
    particles.m_restDensity = fluidRestDensity;
    particles.m_particleRadius = samplingDistance/2;
    particles.m_particleMass = pow(samplingDistance, 3) * fluidRestDensity;

    for (size_t x = 0; x < (size_t)numSamples[0] - 1; x++) {
        for (size_t y = 0; y < (size_t)numSamples[1] - 1; y++) {
            for (size_t z = 0; z < (size_t)numSamples[2] - 1; z++) {
                // NOTE(DENNIS) - Kommen die letzten Samples hier nicht manchmal auch von außerhalb des Volumens?
                Eigen::Vector3d samplingPos ((x + 0.5) * samplingDir[0],
                                             (y + 0.5) * samplingDir[1],
                                             (z + 0.5) * samplingDir[2]);
                samplingPos += bottomLeft;

                particles.positions.push_back(samplingPos);
                particles.velocities.push_back(Eigen::Vector3d::Zero());
                particles.accelerations.push_back(Eigen::Vector3d::Zero());
            }
        }
    }

    return particles;
}
