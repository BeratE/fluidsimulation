#pragma once
#include "particlesystem.h"
#include "boundarysystem.h"
#include "fluidsystem.h"
#include <Eigen/Dense>

namespace learnSPH::System {
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

        FluidSystem sampleFluidBox(Eigen::Vector3d bottomLeft,
                                   Eigen::Vector3d topRight,
                                   double samplingDistance,
                                   double restDensity = 1000.0);

        FluidSystem sampleFluidSDF(std::function<double(Eigen::Vector3d)> const &sdf,
                                   double samplingDistance,
                                   double restDensity = 1000.0);

        FluidSystem sampleFluidSDF(
            std::function<double(Eigen::Vector3d)> const &sdf,            
            Eigen::Vector3d sampleRegionStart,
            Eigen::Vector3d sampleRegionSize,
            double samplingDistance,
            double restDensity = 1000.0);

        
        BoundarySystem sampleBoundaryHollowBox(Eigen::Vector3d bottomLeft,
                                               Eigen::Vector3d topRight,
                                               double samplingDistance,
                                               double restDensity = 1000.0);

        
        BoundarySystem sampleBoundaryPlane(Eigen::Vector3d bottomLeft,
                                           Eigen::Vector3d bottomRight,
                                           Eigen::Vector3d topLeft,
                                           Eigen::Vector3d topRight,
                                           double samplingDistance,
                                           double restDensity = 1000.0);

        BoundarySystem sampleBoundaryMesh(std::string filepath,
                                          double samplingDistance,
                                          double restDensity = 1000.0);
        
        std::vector<Eigen::Vector3d> samplePosTriangle(Eigen::Vector3d a,
                                                    Eigen::Vector3d b,
                                                    Eigen::Vector3d c,
                                                    double samplingDistance);

        std::vector<Eigen::Vector3d> samplePosHollowBox(Eigen::Vector3d bottomLeft,
                                                        Eigen::Vector3d topRight,
                                                        double samplingDistance);

        std::vector<Eigen::Vector3d> samplePosMesh(std::string filepath,
                                                   double samplingDistance);
    };

} // namespace learnSPH
#define Emitter() ParticleEmitter::getInstance()
