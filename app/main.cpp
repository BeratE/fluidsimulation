#include <stdlib.h>     // rand
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>    // std::max

#include <Eigen/Dense>

#include "../learnSPH/util/vtk_writer.h"

int main()
{
    std::cout << "Welcome to the learnSPH framework!!" << std::endl;
    std::cout << "Generating a sample scene...";

    // Generate particles
    std::vector<Eigen::Vector3d> particles;
    for (int i = 0; i < 20; i++) {
        for (int j = 0; j < 20; j++) {
            for (int k = 0; k < 20; k++) {
                double x = static_cast<double>(rand()) / static_cast<float>(RAND_MAX);
                double y = static_cast<double>(rand()) / static_cast<float>(RAND_MAX);
                double z = static_cast<double>(rand()) / static_cast<float>(RAND_MAX);
                particles.push_back(Eigen::Vector3d(x, y, z));
            }
        }
    }

    // Initialize data vectors
    std::vector<double> particles_scalar_data(particles.size());
    std::vector<Eigen::Vector3d> particles_vector_data(particles.size());

    // Simulation loop
    for (int time_step = 0; time_step < 100; time_step++) {

        for (int particle_i = 0; particle_i < (int)particles.size(); particle_i++) {
            // Move particles a bit down in the Z direction
            particles[particle_i][2] -= 0.01;

            // Clamp the Z coord to the floor
            particles[particle_i][2] = std::max(particles[particle_i][2], 0.0);

            // Update scalar data
            particles_scalar_data[particle_i] = particles[particle_i].norm();

            // Update vector data
            particles_vector_data[particle_i] =
                (Eigen::Vector3d(0.5, 0.5, 1.0) - particles[particle_i]).normalized();
        }

        // Save output
        const std::string filename =
            "../res/example_" + std::to_string(time_step) + ".vtk";
        learnSPH::save_particles_to_vtk(filename, particles, particles_scalar_data,
                                        particles_vector_data);
    }

    std::cout << "completed!" << std::endl;
    std::cout << "The scene files have been saved in the folder "
        "`<build_folder>/res`. You can visualize them with Paraview."
              << std::endl;

    return 0;
}
