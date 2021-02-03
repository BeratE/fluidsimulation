#pragma once
#include "solver.h"
#include "system/fluidsystem.h"
#include "system/boundarysystem.h"
#include <vector>
#include <memory>
#include <Eigen/Dense>
#include <CompactNSearch/CompactNSearch.h>
#include <surface/surface.h>

namespace learnSPH {
    class SolverSPH : public Solver {
    public:
        SolverSPH(System::FluidSystem system);
        ~SolverSPH();

        double integrationStep(const std::vector<Eigen::Vector3d> &previousPos) override;
        
        void updateAccelerations(double deltaT);
        
        void updateAccFluidContribution(const size_t i,
                                        const size_t j,
                                        const double ratio_i,
                                        const double ratio_j);

        void updateAccBoundaryContribution(const size_t i,
                                           const size_t k,
                                           const double ratio_i,
                                           System::BoundarySystem& boundary);

    
        // Setter & Getter        
        void setParamStiffness(double val) { m_stiffness = val; }

    private:
        double m_stiffness = 1000.0;
    };
} // namespace learnSPH
