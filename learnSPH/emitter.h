#pragma once
#include "fluidsystem.h"
#include <Eigen/Dense>

namespace learnSPH {

    class Emitter {
    public:
        Emitter() {}
        ~Emitter() {}

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
                                   double fluidRestDensity = 1);
        
    };

} // namespace learnSPH
