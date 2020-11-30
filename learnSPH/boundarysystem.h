#pragma once
#include "particlesystem.h"

namespace learnSPH {
    class BoundarySystem : public ParticleSystem {
        friend class ParticleEmitter;
        
    private:
        BoundarySystem() {}

    public:
        
        void correctVolume();

        // Setter & Getter
        const std::vector<double>& getVolumes() const {return m_volumes; }
        double getParticleVolume(size_t i) const { return m_volumes[i]; }

      private:
        std::vector<double> m_volumes; 
    };
}
