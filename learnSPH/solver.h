#pragma once
#include "fluidsystem.h"
#include "boundarysystem.h"
#include <vector>
#include <CompactNSearch/CompactNSearch.h>
#include <Eigen/Dense>

namespace learnSPH {
    class SolverSPH {
    public:
        SolverSPH(FluidSystem system);
        ~SolverSPH();

        double timeStepCFL();
        double integrationStep();
        void run(std::string file, double milliseconds);
        
        void addBoundary(BoundarySystem boundary);

        // Setter & Getter
        const FluidSystem &getSystem() {return m_system; }

        void setSmoothingEpsilon(double value) {m_smoothEps = value;}

        void setSnapShotAfterMS(double ms) {m_snapShotMS = ms;}
        void setParameterStiffness(double value) {m_system.setStiffness(value); }
        void setParameterViscosity(double value) {m_system.setViscosity(value); }
        void enableGravity(bool value) { m_gravityEnable = value; }
        void enableSmoothing(bool value) {m_smoothingEnable = value; }

      private:
        void semiImplicitEulerStep(double deltaT);
        
        double m_snapShotMS = 20;

        double m_smoothEps = 0.5;
        double m_minTimeStep = 0.002;
        bool m_gravityEnable = true;
        bool m_smoothingEnable = true;
        
        FluidSystem m_system;
        std::vector<BoundarySystem> m_boundaries;
        CompactNSearch::NeighborhoodSearch m_nsearch;

        std::vector<double> m_weightTable;
        std::vector<Eigen::Vector3d> m_weightGradTable;
    };

} // namespace learnSPH
