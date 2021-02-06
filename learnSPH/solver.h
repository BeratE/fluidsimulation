#pragma once
#include "system/fluidsystem.h"
#include "system/boundarysystem.h"
#include <vector>
#include <memory>
#include <Eigen/Dense>
#include <CompactNSearch/CompactNSearch.h>
#include <surface/surface.h>


#define VEC_GRAVITY Eigen::Vector3d(0.0, -9.80665, 0.0)

namespace learnSPH {
class Solver {
public:
    Solver(System::FluidSystem system);
    ~Solver();    
    
    void run(std::string file, double milliseconds);        

    double timeStepCFL();

    size_t addBoundary(const System::BoundarySystem &boundary);
    
    // Setter & Getter
    const System::FluidSystem &getSystem() const { return m_system; }
    const std::vector<System::BoundarySystem> &getBoundaries() const { return m_boundaries; }
    
    double getMaxTimeStepSeconds() { return m_maxTimeStep_s; }
    double getParamSmoothing() { return m_xsphSmoothing; }
    double getSnapShotAfterMS() { return m_snapShotMS; }   
    bool gravityEnabled() { return m_gravityEnable; }
    bool smoothingEnabled() { return m_smoothingEnable; }
    bool tensionEnabled() { return m_tensionEnable; }
    bool adhesionEnabled() { return m_adhesionEnable; }

    void setFluidViscosity(double val) { m_system.setViscosity(val); }
    void setFluidTension(double val) { m_system.setGamma(val); }
    void setBoundaryViscosity(size_t i, double val) {m_boundaries[i].setViscosity(val);}
    void setBoundaryAdhesion(size_t i, double val) { m_boundaries[i].setBeta(val); }
    void setMaxTimeStepSeconds(double val) { m_maxTimeStep_s = val; }
    void setParamSmoothing(double val) { m_xsphSmoothing = val; }
    void setParamDrag(double val) {m_dragCoeff = val; }
    void setSnapShotAfterMS(double ms) { m_snapShotMS = ms; }   
    void setZSortIntervall(unsigned int zSortIntervall) { m_zSortIntervall = zSortIntervall; }
    void enableGravity(bool val) { m_gravityEnable = val; }
    void enableSmoothing(bool val) { m_smoothingEnable = val; }
    void enableTension(bool val) { m_tensionEnable = val; }
    void enableAdhesion(bool val) { m_adhesionEnable = val; }
    void enableDrag(bool val) { m_dragEnable = val; }
    
protected:
    // Specific solvers will overwrite this function. Overhead is minimal.
    virtual void integrationStep(double deltaT) {}
    void semiImplicitEulerStep(double deltaT);
    void initAccelerations();

    void zSort();
    
    double m_snapShotMS = 20;
    double m_maxTimeStep_s = 0.002;
    double m_xsphSmoothing = 0.5;
    double m_dragCoeff = 0.2;
    bool m_gravityEnable = true;
    bool m_smoothingEnable = true;
    bool m_tensionEnable = true;
    bool m_adhesionEnable = true;
    bool m_dragEnable = true;
        
    unsigned int m_zSortIntervall = 100;
    System::FluidSystem m_system;
    std::vector<System::BoundarySystem> m_boundaries;
    std::shared_ptr<CompactNSearch::NeighborhoodSearch> mp_nsearch; 
};
} // namespace learnSPH
