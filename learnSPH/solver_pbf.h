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
    class SolverPBF : public Solver {
    public:
        SolverPBF(System::FluidSystem system);
        ~SolverPBF();

        double integrationStep(const std::vector<Eigen::Vector3d> &previousPos) override;
        
        void updateAccelerations(const double deltaT);
    
        void updateAccFluidContribution(const size_t i,
                                        const size_t j,
                                        const double ratio_i,
                                        const double ratio_j);

        void updateAccBoundaryContribution(const size_t i,
                                           const size_t k,
                                           const double ratio_i,
                                           System::BoundarySystem& boundary);

        void updatePositionsWithConstraints();
    
        // Setter & Getter
        void setNumIterations(size_t value) {m_npbfIterations = value;}

    private:
        size_t m_npbfIterations = 3;
        double m_eps = 0.0001;
    };
} // namespace learnSPH
