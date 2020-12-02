#include "boundarysystem.h"

using namespace learnSPH;
using namespace learnSPH::Kernel;
using namespace learnSPH::System;
using namespace CompactNSearch;

BoundarySystem::BoundarySystem(size_t size, bool fill)
    : ParticleSystem(size, fill)
{
    m_volumes.resize(size);
    if (fill) {
        std::fill(m_volumes.begin(), m_volumes.end(), 0.0);
    }
}

double BoundarySystem::particleMass(size_t i) const
{
    return m_restDensity * m_volumes[i];
}

void BoundarySystem::updateVolume()
{   
    NeighborhoodSearch nsearch(CubicSpline::support(smoothingLength()));
    addToNeighborhood(nsearch);
    nsearch.find_neighbors();

    PointSet const& ps = nsearch.point_set(getPointSetID());
    m_volumes.resize(ps.n_points());
    
    for (size_t i = 0; i < ps.n_points(); i++) {
        double vol = CubicSpline::weight(m_positions[i], m_positions[i],
                                         smoothingLength());
        for (size_t j = 0; j < ps.n_neighbors(getPointSetID(), i); j++) {
            const size_t nid = ps.neighbor(getPointSetID(), i, j);
            vol += CubicSpline::weight(m_positions[i], m_positions[nid],
                                       smoothingLength());
        }
        m_volumes[i] = 1.0 / vol;
    }
}
