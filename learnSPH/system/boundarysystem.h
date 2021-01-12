#pragma once
#include "particlesystem.h"

namespace learnSPH::System {
    class BoundarySystem : public ParticleSystem {
        friend class ParticleEmitter;
        
    public:
        BoundarySystem(double radius, double density, size_t size, bool fill = true);

        void updateVolumes();
        void transform(Eigen::Matrix4d transform);

        
        double getParticleMass(size_t i) const;       

        // Setter & Getter
        double getParticleVolume(size_t i) const { return m_volumes[i]; }
        
        const std::vector<double>& getVolumes() const {return m_volumes; }

      private:
        std::vector<double> m_volumes; 
    };
} // namespace learnSPH::System
