#include "boundarysystem.h"

using namespace learnSPH;
using namespace learnSPH::Kernel;
using namespace learnSPH::System;
using namespace CompactNSearch;

BoundarySystem::BoundarySystem(double radius, double density, size_t size, bool fill)
    : ParticleSystem(radius, density, size, fill)
{
    m_volumes.resize(size);
    if (fill) {
        std::fill(m_volumes.begin(), m_volumes.end(), 0.0);
    }
}

double BoundarySystem::getParticleMass(size_t i) const
{
    return m_restDensity * m_volumes[i];
}

void BoundarySystem::updateVolumes()
{
    NeighborhoodSearch nsearch(CubicSpline::support(m_smoothingLength));
    size_t id = nsearch.add_point_set(m_positions.front().data(), m_positions.size());
    nsearch.find_neighbors();

    PointSet const& ps = nsearch.point_set(id);
    m_volumes.resize(ps.n_points());
    
    for (size_t i = 0; i < ps.n_points(); i++) {
        double vol = CubicSpline::weight(m_positions[i], m_positions[i],
                                         m_smoothingLength);
        for (size_t j = 0; j < ps.n_neighbors(id, i); j++) {
            const size_t nid = ps.neighbor(id, i, j);
            vol += CubicSpline::weight(m_positions[i], m_positions[nid],
                                       m_smoothingLength);
        }
        m_volumes[i] = 1.0 / vol;
    }
}

void BoundarySystem::transform(Eigen::Matrix4d transform)
{
    for (auto &p : m_positions) {
        Eigen::Vector4d tmp = Eigen::Vector4d(p(0), p(1), p(2), 1.0);
        p = (transform * tmp).segment(0, 3);
    }
}
