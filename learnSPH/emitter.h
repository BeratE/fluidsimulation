#pragma once
#include "particlesystem.h"
#include "boundarysystem.h"
#include "fluidsystem.h"
#include <Eigen/Dense>

namespace learnSPH {
    /* The Emitter (singleton) is responsible for generating and manipulating 
       particle systems. */
    class ParticleEmitter {        
    private:
        ParticleEmitter() {}
        ParticleEmitter(const ParticleEmitter&) = delete;
        ParticleEmitter &operator=(ParticleEmitter) = delete;

    public:
        static ParticleEmitter &getInstance()
        {
            static ParticleEmitter instance;
            return instance;
        }

        /* Uniformly sample an axis box with fluid particles.
         * @param bottomLeft       - bottom left corner of 3d axis aligned box eg.
         * (0,0,0).
         * @param topRight         - top right corner diagonal to bottom left eg.
         * (1,1,1).
         * @param fluidRestDensity - rest density of particle system.
         * @param samplingDistance - distance between particle samples.
         *                           Effective particle diameter.
         * @return FuidSystem particles in the axis aligned box. */
        FluidSystem sampleFluidBox(Eigen::Vector3d bottomLeft,
                                   Eigen::Vector3d topRight,
                                   double samplingDistance,
                                   double restDensity = 1000.0);

        
        /* Construct a boundary from the given positions, rest density and mass */
        BoundarySystem sampleBoundaryHollowBox(Eigen::Vector3d bottomLeft,
                                               Eigen::Vector3d topRight,
                                               double samplingDistance,
                                               double restDensity = 1000.0);

        
        /* Sample positions in a triangle specified by 3d points a, b , c with
         * sampling distance apart from each other (hexagonal).
         * @param a, b , c         - 3d points of the triangle.
         * @param samplingDistance - distance the sample points are apart.
         * @return Vector of sampled 3d points within the triangle. */
        std::vector<Eigen::Vector3d> samplePosTriangle(Eigen::Vector3d a,
                                                    Eigen::Vector3d b,
                                                    Eigen::Vector3d c,
                                                    double samplingDistance);

        /* Sample a hollow box spanning from the bottom left to the top right corner.
         * @param bottomLeft - bottom left box corner.
         * @param topRight   - top right box corner diagonal to bottom left.
         * @param samplingDistance - particles are sampling distance apart.
         * @return vector of positions of sample points. */
        std::vector<Eigen::Vector3d> samplePosHollowBox(Eigen::Vector3d bottomLeft,
                                                        Eigen::Vector3d topRight,
                                                        double samplingDistance);
    };

} // namespace learnSPH
#define Emitter() ParticleEmitter::getInstance()
