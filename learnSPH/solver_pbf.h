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
    
    double integrationStep();
    void run(std::string file, double milliseconds,
             std::vector<Surface::SurfaceInformation>* pOutSurfaceInfos = nullptr) override;

    void newRun(std::string file, double milliseconds,
        std::vector<Surface::SurfaceInformation>* pOutSurfaceInfos = nullptr);
    double newIntegrationStep();
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
    double C(size_t i);
    double S(size_t i);
    Eigen::Vector3d deltaX(size_t i, std::vector<double> lambda);

    size_t m_npbfIterations = 3;
    double m_eps = 0.0001;
};
} // namespace learnSPH
