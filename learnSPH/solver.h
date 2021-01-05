#pragma once
#include "system/fluidsystem.h"
#include "system/boundarysystem.h"
#include <vector>
#include <memory>
#include <Eigen/Dense>
#include <CompactNSearch/CompactNSearch.h>
#include <surface/surface.h>


#define VEC_GRAVITY Eigen::Vector3d(0.0, -0.980665, 0.0)

namespace learnSPH {
class Solver {
public:
    Solver(System::FluidSystem system);
    ~Solver();

    double timeStepCFL();
    //double integrationStep() {}
    void run(std::string file, double milliseconds) {}

    size_t addBoundary(const System::BoundarySystem &boundary);

    // Setter & Getter
    const System::FluidSystem &getSystem() const { return m_system; }

    void setSnapShotAfterMS(double ms) { m_snapShotMS = ms; }
    void setFluidViscosity(double val) { m_system.setViscosity(val); }
    void setBoundaryViscosity(size_t i, double val) {m_boundaries[i].setViscosity(val);}
    void enableGravity(bool val) { m_gravityEnable = val; }
    void enableSmoothing(bool val) { m_smoothingEnable = val; }

protected:        
    void applyExternalForces();
    //void semiImplicitEulerStep(double deltaT) {}

    double m_snapShotMS = 20;
    double m_maxTimeStep_s = 0.002;
    bool m_gravityEnable = true;
    bool m_smoothingEnable = true;
        
    System::FluidSystem m_system;
    std::vector<System::BoundarySystem> m_boundaries;
    std::shared_ptr<CompactNSearch::NeighborhoodSearch> mp_nsearch;
};
} // namespace learnSPH
