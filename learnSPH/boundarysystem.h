#pragma once
#include "particlesystem.h"

namespace learnSPH {
    class BoundarySystem : public ParticleSystem {
        friend class Emitter;
        
    private:
      BoundarySystem() {}

    public:
        
        void correctVolume();
        double getParticleVolume(size_t i) { return m_volumes[i]; }

      private:
        std::vector<double> m_volumes; 
    };
}
