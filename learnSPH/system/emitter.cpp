#include "emitter.h"
#include <math.h>

using namespace learnSPH;
using namespace learnSPH::System;
using namespace Eigen;

FluidSystem ParticleEmitter::sampleFluidBox(Vector3d bottomLeft,
                                            Vector3d topRight,
                                            double samplingDistance,
                                            double restDensity)
{
    const Eigen::Vector3d diagonal = topRight - bottomLeft;
    const Eigen::Vector3d diagonalSign = diagonal.array() / diagonal.array().abs();
    const Eigen::Vector3d samplingDir = samplingDistance * diagonalSign;
    const Eigen::Vector3d numSamples = Eigen::floor(Eigen::abs(diagonal.array())
                                                    / samplingDistance);
        
    size_t size = numSamples[0] * numSamples[1] * numSamples[2];
    
    FluidSystem particles(samplingDistance/2, restDensity, size);
    
    for (size_t x = 0; x < (size_t)numSamples[0]; x++) {
        for (size_t y = 0; y < (size_t)numSamples[1]; y++) {
            for (size_t z = 0; z < (size_t)numSamples[2]; z++) {
                Eigen::Vector3d samplingPos  = bottomLeft
                    + Eigen::Vector3d((x + 0.5) * samplingDir[0],
                                      (y + 0.5) * samplingDir[1],
                                      (z + 0.5) * samplingDir[2]);
                size_t index =
                    x*numSamples[1]*numSamples[2] +
                    y*numSamples[2] +
                    z;

                particles.m_positions[index] = samplingPos;
            }
        }
    }

    return particles;
}

BoundarySystem ParticleEmitter::sampleBoundaryHollowBox(Vector3d bottomLeft,
                                                        Vector3d topRight,
                                                        double samplingDistance,
                                                        double restDensity)
{
    auto pos = samplePosHollowBox(bottomLeft, topRight, samplingDistance*0.8);
    
    BoundarySystem boundary(samplingDistance/2, restDensity, pos.size());
    boundary.m_positions = pos;
    boundary.updateVolumes();
    return boundary;
}

BoundarySystem ParticleEmitter::sampleBoundaryPlane(Eigen::Vector3d bottomLeft,
                                                    Eigen::Vector3d bottomRight,
                                                    Eigen::Vector3d topLeft,
                                                    Eigen::Vector3d topRight,
                                                    double samplingDistance,
                                                    double restDensity)
{
    auto triangle1 = samplePosTriangle(bottomLeft, bottomRight, topRight, samplingDistance*0.8);
    auto triangle2 = samplePosTriangle(topRight, topLeft, bottomLeft, samplingDistance*0.8);
    triangle1.insert(triangle1.end(), triangle2.begin(), triangle2.end());
    
    BoundarySystem boundary(samplingDistance/2, restDensity, triangle1.size());
    boundary.m_positions = triangle1;
    boundary.updateVolumes();
    return boundary;
}

std::vector<Vector3d> ParticleEmitter::samplePosTriangle(Vector3d a,
                                                         Vector3d b,
                                                         Vector3d c,
                                                         double samplingDistance)
{   
    // Construct orthonormal base of the 2d plane with origin at a spanned by the triangle
    Vector3d u = (b-a).normalized();
    Vector3d n = u.cross((c-a).normalized()).normalized();
    Vector3d v = n.cross(u).normalized();
    Matrix<double, 3, 4> planeBasis;
    planeBasis.row(0) = Vector4d(u(0), u(1), u(2), -u.dot(a));
    planeBasis.row(1) = Vector4d(v(0), v(1), v(2), -v.dot(a));
    planeBasis.row(2) = Vector4d(0.0, 0.0, 0.0, 1.0);
    
    // Convert triangle points to new 2d basis (discard homogenous coordinate)
    Vector2d x = (planeBasis * a.homogeneous()).segment(0, 2);
    Vector2d y = (planeBasis * b.homogeneous()).segment(0, 2);
    Vector2d z = (planeBasis * c.homogeneous()).segment(0, 2);

    // Return the hessian line function of the edge from point a to point b,
    // offset by half the sampling distance from the opposite triangle point o.
    auto constructEdge = [samplingDistance](Eigen::Vector2d a,
                                            Eigen::Vector2d b,
                                            Eigen::Vector2d o)
        -> std::function<double(Eigen::Vector2d)> {
        Eigen::Vector2d normal = (b-a);
        double tmp = normal(0);
        normal(0) = normal(1);
        normal(1) = -tmp;
        normal.normalize();
        Eigen::Vector2d offset = ((a+(b-a)/2)-o).normalized()*samplingDistance/2;
        double distance = normal.dot(a + offset);
        return [=](Eigen::Vector2d point) mutable {
            return normal.dot(point) - distance;
        };
    };

    auto edgeXY = constructEdge(x, y, z);
    auto edgeYZ = constructEdge(y, z, x);
    auto edgeZX = constructEdge(z, x, y);
    
    // Generate points in the 2d rectancle enclosing the triangle with the
    // pineda algorithm
    Eigen::Vector2d rectMin = -(samplingDistance/2.0) *Vector2d(1.0, 1.0); 
    Eigen::Vector2d rectMax =  (samplingDistance/2.0)*Vector2d(1.0, 1.0); 
    rectMin(0) += (x(0) < y(0)) ? (x(0) < z(0) ? x(0) : z(0)) : (y(0) < z(0) ? y(0) : z(0));
    rectMin(1) += (x(1) < y(1)) ? (x(1) < z(1) ? x(1) : z(1)) : (y(1) < z(1) ? y(1) : z(1));
    rectMax(0) += (x(0) > y(0)) ? (x(0) > z(0) ? x(0) : z(0)) : (y(0) > z(0) ? y(0) : z(0));
    rectMax(1) += (x(1) > y(1)) ? (x(1) > z(1) ? x(1) : z(1)) : (y(1) > z(1) ? y(1) : z(1));

    // Sample points in a hexagonal raster of
    // integer steps of sampling distance in j direction and
    // sampling distance/sqrt(2) in i direction
    std::vector<Vector2d> samples;
    const size_t nStepsI = std::floor((rectMax(0)-rectMin(0))/(samplingDistance/sqrt(2)));
    const size_t nStepsJ = std::floor((rectMax(1)-rectMin(1))/samplingDistance);
    for (size_t i = 0; i < nStepsI; i++) {
        for (size_t j = 0; j < nStepsJ; j++) {
            Vector2d sample(rectMin(0) + i * (samplingDistance/sqrt(2)),
                            rectMin(1) + j * samplingDistance);
            
            sample(1) += ((i%2) * samplingDistance*sqrt(2));
            
            if (edgeXY(sample) <= 0 && edgeYZ(sample) <= 0 && edgeZX(sample) <= 0)
                samples.push_back(sample);
        }
    }

    std::vector<Vector3d> positions;
    // Transform sample points back to 3d standard basis
    for (const auto &s : samples) {
        Vector3d sample_point = (planeBasis.transpose()
                                 * s.homogeneous()).segment(0, 3) + a;
        positions.push_back(sample_point);
    }

    return positions;
}


std::vector<Vector3d> ParticleEmitter::samplePosHollowBox(Vector3d bottomLeft,
                                                          Vector3d topRight,
                                                          double samplingDistance)
{
    std::vector<Vector3d> positions;
    // Generate corner vertices
    Vector3d dir = topRight - bottomLeft;
    std::vector<Vector3d> vertices;
    vertices.push_back(Vector3d(bottomLeft(0), bottomLeft(1), bottomLeft(2)));
    vertices.push_back(Vector3d(bottomLeft(0), bottomLeft(1)+dir(1), bottomLeft(2)));
    vertices.push_back(Vector3d(bottomLeft(0)+dir(0), bottomLeft(1), bottomLeft(2)));
    vertices.push_back(Vector3d(bottomLeft(0)+dir(0), bottomLeft(1)+dir(1), bottomLeft(2)));
    vertices.push_back(Vector3d(bottomLeft(0), bottomLeft(1), bottomLeft(2)+dir(1)));
    vertices.push_back(Vector3d(bottomLeft(0), bottomLeft(1)+dir(1), bottomLeft(2)+dir(1)));
    vertices.push_back(Vector3d(bottomLeft(0)+dir(0), bottomLeft(1),bottomLeft(2)+dir(1)));
    vertices.push_back(Vector3d(bottomLeft(0)+dir(0), bottomLeft(1)+dir(1), bottomLeft(2)+dir(1)));

    std::vector<std::vector<Vector3d>> v;
    v.push_back(samplePosTriangle(vertices[2], vertices[1], vertices[3], samplingDistance));
    v.push_back(samplePosTriangle(vertices[0], vertices[1], vertices[2], samplingDistance));
    v.push_back(samplePosTriangle(vertices[0], vertices[4], vertices[2], samplingDistance));
    v.push_back(samplePosTriangle(vertices[2], vertices[4], vertices[6], samplingDistance));
    v.push_back(samplePosTriangle(vertices[0], vertices[4], vertices[1], samplingDistance));
    v.push_back(samplePosTriangle(vertices[1], vertices[4], vertices[5], samplingDistance));
    v.push_back(samplePosTriangle(vertices[1], vertices[5], vertices[3], samplingDistance));
    v.push_back(samplePosTriangle(vertices[3], vertices[5], vertices[7], samplingDistance));
    v.push_back(samplePosTriangle(vertices[2], vertices[6], vertices[3], samplingDistance));
    v.push_back(samplePosTriangle(vertices[3], vertices[6], vertices[7], samplingDistance));
    v.push_back(samplePosTriangle(vertices[4], vertices[5], vertices[6], samplingDistance));
    v.push_back(samplePosTriangle(vertices[6], vertices[5], vertices[7], samplingDistance));
    
    for (const auto &p : v) {
        positions.insert(positions.end(), p.begin(), p.end());
    }

    return positions;
}
