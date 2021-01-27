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
    
    void newRun(std::string file, double milliseconds, std::vector<Surface::SurfaceInformation>* pOutSurfaceInfos = nullptr);
    void newSemiImplicitEulerStep(double deltaT);

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

    double integrationStep();
    double newIntegrationStep();

    void run(std::string file, double milliseconds,
             std::vector<Surface::SurfaceInformation>* pOutSurfaceInfos = nullptr) override;
    
    // Setter & Getter        
    void setParamStiffness(double val) { m_stiffness = val; }
    void setParameterDrag(double val) { m_drag = val; }    

  private:
    double m_drag = 0.2;    
    double m_stiffness = 1000.0;
};
} // namespace learnSPH
