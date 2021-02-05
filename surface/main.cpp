#include <iostream>
#include <Eigen/Dense>
#include <vtkio/VTKFile.h>
#include <CompactNSearch/CompactNSearch.h>
#include "vtk_writer.h"

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

    std::string filename(argv[1]);
    std::string comments;
    std::vector<Eigen::Vector3d> positions;
    
    learnSPH::readParticlesFromVTK(filename, positions, comments);

    // Get Properties for lookup table    
    size_t numBins;
    double smoothingLength;
    std::string stra, strb;
    size_t i;
    i = getNumStr(comments, 0, stra);
    if (i == -1) return -1;
    i = getNumStr(comments, i, strb);
    if (i == -1) return -1;
    numBins = std::stoi(stra);
    smoothingLength = std::stod(strb);

    double searchradius = learnSPH::Kernel::CubicSpline::support(smoothingLength);
    CompactNSearch::NeighborhoodSearch nsearch(searchradius);
    learnSPH::Kernel::CubicSpline::Table splineLut(smoothingLength, numBins);

    
    std::vector<double> densities(positions.size());
    for (size_t i = 0; i < densities.size(); i++) {
        
    }

    

    return 0;
}
