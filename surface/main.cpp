#include <iostream>
#include <Eigen/Dense>
#include <vtkio/VTKFile.h>
#include <CompactNSearch/CompactNSearch.h>
#include "vtk_writer.h"
#include "learnSPH/kernel.h"
#include "surface.h"

using namespace learnSPH::Surface;

size_t getNumStr(std::string str, size_t i, std::string &outstr)
{
    size_t a = str.find(':', i);
    size_t b = str.find(';', a+1);

    if (a == std::string::npos || b == std::string::npos) {
        std::cout << "Error, .vtk file does not contain proper comment" << std::endl;
        return -1;
    }

    outstr = str.substr(a+1, b-a-1);
    return b+1;
}

int main(int argc, char *argv[])
{
    if (argc < 2) {
        std::cout << "Error, .vtk input file required" << std::endl;
    }

    const std::string filename(argv[1]);

    std::cout << "Reading file " << filename << std::endl;

    // Read data from VTK file
    std::string comments;
    std::vector<Eigen::Vector3d> positions;    
    learnSPH::readParticlesFromVTK(filename, positions, comments);

    // Reconstruct properties  
    size_t numBins;
    double smoothingLength;
    std::string stra, strb;
    size_t i;
    i = getNumStr(comments, 0, stra);
    if (i == -1) return -1;
    i = getNumStr(comments, i, strb);
    if (i == -1) return -1;
    smoothingLength = std::stod(stra);
    numBins = std::stoi(strb);

    learnSPH::Kernel::CubicSpline::Table splineLut(smoothingLength, numBins);

    // Compute neighborhood
    double searchradius = learnSPH::Kernel::CubicSpline::support(smoothingLength);
    CompactNSearch::NeighborhoodSearch nsearch(searchradius);
    size_t pid = nsearch.add_point_set(positions.front().data(),
                                        positions.size());
    nsearch.find_neighbors();
    CompactNSearch::PointSet const ps = nsearch.point_set(pid);

    // Compute normalised densities
    std::vector<double> densities(positions.size());

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < ps.n_points(); i++) {
        double normDensity = splineLut.weight(positions[i], positions[i]);
        for (size_t j = 0; j < ps.n_neighbors(pid, i); j++) {
            const size_t k = ps.neighbor(pid, i, j);
            normDensity += splineLut.weight(positions[i], positions[k]);
        }
        densities[i] = normDensity;        
    }    
    
    // Reconstruct SDF
    const double paramC = 0.6;
    const double ratioSmoothingLengthSamplingStep = 2.0;
    std::vector<Eigen::Vector3d> gridVerts;
    std::vector<double> gridSDF;    
    Eigen::Vector3i gridDims;
    discretizeFluidSystemSDF(
        positions, densities, splineLut, smoothingLength,
        smoothingLength / ratioSmoothingLengthSamplingStep, paramC,
        &gridSDF, &gridVerts, &gridDims);

    // Extract Surface
    std::vector<Eigen::Vector3d> vertices;
    std::vector<std::array<int, 3>> triangles;
    marchCubes(gridDims, gridSDF, gridVerts, vertices, triangles);

    // Split into path and name
    size_t m = filename.find_last_of('/');
    std::string path = filename.substr(0, m != std::string::npos ? m : 0);
    std::string name = filename.substr(m != std::string::npos ? m : 0, std::string::npos);

    // Find optional filenumber
    std::string filenumber = "";
    char const* digits = "0123456789";    
    std::size_t n = name.find_first_of(digits);
    if (n != std::string::npos) {
        std::size_t const m = name.find_first_not_of(digits, n);
        filenumber = name.substr(n, m != std::string::npos ? m-n : m);
    }    

    // Output file
    std::stringstream fn;
    fn << path << name.substr(0, n) << "_surface" << filenumber << ".vtk";

    learnSPH::writeMeshToVTK(fn.str(), vertices, triangles);
    
    return 0;
}
