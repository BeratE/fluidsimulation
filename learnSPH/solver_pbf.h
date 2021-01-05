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
    void run(std::string file, double milliseconds, std::vector<Surface::SurfaceInformation>* pOutSurfaceInfos = nullptr);
    
    // Setter & Getter        

  private:
    double C(size_t i);
    double S(size_t i);
    Eigen::Vector3d deltaX(size_t i, std::vector<double> lambda);
    void semiImplicitEulerStep(double deltaT);

    size_t m_npbfIterations = 3;
};
} // namespace learnSPH
