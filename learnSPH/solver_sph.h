#pragma once
#include "solver.h"
#include "system/fluidsystem.h"
#include "system/boundarysystem.h"
#include <vector>
#include <memory>
#include <Eigen/Dense>
#include <CompactNSearch/CompactNSearch.h>

namespace learnSPH {
struct ParameterSPH {
    double drag = 0.2;
    double stiffness = 1000.0;
    double smoothing = 0.5;
};

class SolverSPH : public Solver {
  public:
    SolverSPH(System::FluidSystem system);
    ~SolverSPH();
    
    double integrationStep();      
    void run(std::string file, double milliseconds);
    
    // Setter & Getter        
    void setParamSmoothing(double val) { m_param.smoothing = val; }
    void setParamStiffness(double val) { m_param.stiffness = val; }
    void setParameterDrag(double val) { m_param.drag = val; }     
    
  private:        
    void semiImplicitEulerStep(double deltaT);                          
    
    ParameterSPH m_param;
};
} // namespace learnSPH
