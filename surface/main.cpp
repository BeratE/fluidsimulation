#include <iostream>
#include <Eigen/Dense>
#include <vtkio/VTKFile.h>

int main(int argc, char *argv[])
{
    if (argc <= 1) {
        std::cout << "Error, .vtk input file required" << std::endl;
    }

    vtkio::VTKFile vtk_file;
    
    vtk_file.read(argv[0]);
   
    size_t numPoints = vtk_file.get_number_of_points();
    std::vector<Eigen::Vector3d> positions(numPoints);
    std::vector<double> densities(numPoints);

    vtk_file.get_points_to_twice_indexable(positions);
    vtk_file.get_point_data_to_indexable("scalar", densities);
            
    return 0;
}
