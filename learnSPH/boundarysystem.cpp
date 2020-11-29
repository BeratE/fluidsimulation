#include "boundarysystem.h"
#include "kernel.h"

using namespace learnSPH;
using namespace learnSPH::Kernel;
using namespace CompactNSearch;

void BoundarySystem::correctVolume()
{   
    NeighborhoodSearch nsearch(CubicSpline::support(smoothingLength()));

    unsigned int id  = nsearch.add_point_set(m_positions.front().data(),
                                             m_positions.size());
    nsearch.find_neighbors();

    PointSet const& ps = nsearch.point_set(id);
    m_volumes.resize(ps.n_points());
    
    for (size_t i = 0; i < ps.n_points(); i++) {
        double vol = CubicSpline::weight(m_positions[i], m_positions[i],
                                         smoothingLength());
        for (size_t j = 0; j < ps.n_neighbors(id, i); j++) {
            const size_t nid = ps.neighbor(id, i, j);
            vol += CubicSpline::weight(m_positions[i], m_positions[nid],
                                       smoothingLength());
        }
        m_volumes[i] = 1.0 / vol;
    }
}
