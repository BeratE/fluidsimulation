#pragma once
#include "system/fluidsystem.h"
#include "system/boundarysystem.h"
#include <vector>
#include <memory>
#include <Eigen/Dense>
#include <CompactNSearch/CompactNSearch.h>


#define VEC_GRAVITY Eigen::Vector3d(0.0, -0.980665, 0.0)

namespace learnSPH {
    
    struct ParameterSPH {
        double drag = 0.2;
        double stiffness = 1000.0;
        double smoothing = 0.5;
    };

    class SolverSPH {
    public:
        SolverSPH(System::FluidSystem system);
        ~SolverSPH();

        double timeStepCFL();
        double integrationStep();
        
        void run(std::string file, double milliseconds);

        size_t addBoundary(System::BoundarySystem boundary);

        // Setter & Getter
        const System::FluidSystem &getSystem() const { return m_system; }

        void setSnapShotAfterMS(double ms) { m_snapShotMS = ms; }
        
        void setParamSmoothing(double val) { m_param.smoothing = val; }
        void setParamStiffness(double val) { m_param.stiffness = val; }
        void setParameterDrag(double val) { m_param.drag = val; }

        void setFluidViscosity(double val) {m_system.setViscosity(val);}
        void setBoundaryViscosity(size_t i, double val) {m_boundaries[i].setViscosity(val);}

        void enableGravity(bool val) { m_gravityEnable = val; }
        void enableSmoothing(bool val) { m_smoothingEnable = val; }

    private:
        void applyExternalForces();
        void semiImplicitEulerStep(double deltaT);

        double m_snapShotMS = 20;

        ParameterSPH m_param;
        
        double m_maxTimeStep = 0.002; // millisecods
        bool m_gravityEnable = true;
        bool m_smoothingEnable = true;

        System::FluidSystem m_system;
        std::vector<System::BoundarySystem> m_boundaries;
        std::shared_ptr<CompactNSearch::NeighborhoodSearch> mp_nsearch;
    };

} // namespace learnSPH
