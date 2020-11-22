#include <stdlib.h>     // rand
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>    // std::max
#include <Eigen/Dense>
#include "../learnSPH/kernel.h"
#include "../learnSPH/particlesystem.h"
#include "../learnSPH/util/vtk_writer.h"

using namespace learnSPH;

int main()
{
    std::cout << "Welcome to the learnSPH framework!!" << std::endl;
    std::cout << "Generating a sample scene...";

    // Generate particles
    const double restDensity = 1.0;
    const double particleDiameter = 0.1;
    const double smoothingLength = particleDiameter * Kernel::Parameter::TUNING;
    // Sample Particles in a Box
    ParticleSystem::FluidSystem particles
        = ParticleSystem::sampleBox(Eigen::Vector3d(0, 0, 0),
                                    Eigen::Vector3d(1, 1, 1),
                                    restDensity, particleDiameter);
    // Compute neighborhood information of fluid particles
    CompactNSearch::NeighborhoodSearch nsearch(Kernel::CubicSpline::support(smoothingLength));
    particles.id = nsearch.add_point_set(particles.positions.front().data(),
                                         particles.positions.size());
    nsearch.find_neighbors();
    
    ParticleSystem::estimateDensity(particles, nsearch);

    const std::string filename = "../../res/test_p1_d005_n12.vtk";
    learnSPH::save_particles_to_vtk(filename, particles.positions, particles.densities);

    
    // // Initialize data vectors
    // std::vector<double> particles_scalar_data(particles.positions.size());
    // std::vector<Eigen::Vector3d> particles_vector_data(particles.positions.size());

    // // Simulation loop
    // for (int time_step = 0; time_step < 100; time_step++) {

    //     for (int particle_i = 0; particle_i < (int)particles.positions.size(); particle_i++) {
    //         // Move particles a bit down in the Z direction
    //         particles.positions[particle_i][2] -= 0.01;

    //         // Clamp the Z coord to the floor
    //         particles.positions[particle_i][2] = std::max(particles.positions[particle_i][2], 0.0);

    //         // Update scalar data
    //         particles_scalar_data[particle_i] = particles.positions[particle_i].norm();

    //         // Update vector data
    //         particles_vector_data[particle_i] =
    //             (Eigen::Vector3d(0.5, 0.5, 1.0) - particles.positions[particle_i]).normalized();
    //     }

    //     // Save output
    //     const std::string filename =
    //         "../res/example_" + std::to_string(time_step) + ".vtk";
    //     learnSPH::save_particles_to_vtk(filename, particles.positions, particles_scalar_data,
    //                                     particles_vector_data);
    // }

    
    std::cout << "completed!" << std::endl;
    std::cout << "The scene files have been saved in the folder "
        "`<build_folder>/res`. You can visualize them with Paraview."
              << std::endl;

    return 0;
}
