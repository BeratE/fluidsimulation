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
        double getBeta() const { return m_beta; }
        const std::vector<double>& getVolumes() const {return m_volumes; }
        
        void setBeta(const double beta) { m_beta = beta; }

      private:
        std::vector<double> m_volumes; 
        double m_beta = 1.0;
    };
} // namespace learnSPH::System
