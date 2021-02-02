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

    void newRun(
        std::string file, double milliseconds,
        std::vector<Surface::SurfaceInformation>* pOutSurfaceInfos = nullptr) override;
    
    double newIntegrationStep(std::vector<Eigen::Vector3d>& previousPos);
    void newSemiImplicitEulerStep(const double deltaT);
    void updateAccFluidContribution(std::vector<Eigen::Vector3d>& accelerations,
        const size_t i,
        const size_t j,
        const double ratio_i,
        const double ratio_j);

    void updateAccBoundaryContribution(std::vector<Eigen::Vector3d>& accelerations,
        const size_t i,
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
