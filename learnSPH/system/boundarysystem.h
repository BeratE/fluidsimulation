#pragma once
#include "particlesystem.h"

namespace learnSPH::System {
    class BoundarySystem : public ParticleSystem {
        friend class ParticleEmitter;
        
    public:
        BoundarySystem() {}
        BoundarySystem(size_t size, bool fill = true);

        double particleMass(size_t i) const;
        
        void updateVolume();

        // Setter & Getter
        double getParticleVolume(size_t i) const { return m_volumes[i]; }
        
        const std::vector<double>& getVolumes() const {return m_volumes; }

      private:
        std::vector<double> m_volumes; 
    };
} // namespace learnSPH::System
