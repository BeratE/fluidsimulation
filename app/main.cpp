#include <stdlib.h>     // rand
#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <cmath>
#include <cfloat>
#include <algorithm>    // std::max

#include <CompactNSearch/CompactNSearch>
#include <Eigen/Dense>

#include "../learnSPH/kernel.h"
#include "../learnSPH/util/vtk_writer.h"

std::vector<std::vector<unsigned int>> bruteForceNeighborhoodSearch(const std::vector<std::array<double, 3>> positions, const double radius) {
    std::size_t positions_size = positions.size();
    std::vector<std::vector<unsigned int>> neighborhoods(positions_size);

    for (std::size_t i = 0; i < positions_size; i++) {
        for (std::size_t j = 0; j < positions_size; j++) {
            if (i == j) {
                continue;
            }
            Eigen::Vector3d pos_i(positions[i].data());
            Eigen::Vector3d pos_j(positions[j].data());

            if ((pos_i - pos_j).norm() < radius) {
                neighborhoods[i].push_back(j);
            }

        }
    }

    return neighborhoods;
}

std::vector<std::array<double, 3>> generateParticles(const unsigned int particlesPerDimension, const double sizeOfDimension) {
    std::vector<std::array<double, 3>> particles;
    for (int i = 0; i < particlesPerDimension; i++) {
        for (int j = 0; j < particlesPerDimension; j++) {
            for (int k = 0; k < particlesPerDimension; k++) {
                double x = static_cast<double> (sizeOfDimension * (double)rand()) / static_cast <float> (RAND_MAX);
                double y = static_cast<double> (sizeOfDimension * (double)rand()) / static_cast <float> (RAND_MAX);
                double z = static_cast<double> (sizeOfDimension * (double)rand()) / static_cast <float> (RAND_MAX);
                particles.push_back(std::array<double, 3>{x, y, z});
            }
        }
    }
    return particles;
}



void compareNeighborhoodSearches(const unsigned int particlesPerDimension, const double sizeOfDimension, const double radius) {

    for (unsigned int curParticlesPerDimension = 1; curParticlesPerDimension <= particlesPerDimension; curParticlesPerDimension++) {
        std::vector<std::array<double, 3>> particles = generateParticles(curParticlesPerDimension, sizeOfDimension);

        auto t1 = std::chrono::high_resolution_clock::now();
        std::vector<std::vector<unsigned int>> bruteForceNeighborhood = bruteForceNeighborhoodSearch(particles, radius);
        auto t2 = std::chrono::high_resolution_clock::now();
        auto durationBruteForce = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

        CompactNSearch::NeighborhoodSearch nsearch(radius);
        unsigned int pointSetId = nsearch.add_point_set(particles.front().data(), particles.size());
        t1 = std::chrono::high_resolution_clock::now();
        nsearch.find_neighbors();
        t2 = std::chrono::high_resolution_clock::now();
        auto durationCompactNSearch = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

        CompactNSearch::PointSet const& ps = nsearch.point_set(pointSetId);

        for (std::size_t i = 0; i < particles.size(); i++) {
            assert(bruteForceNeighborhood[i].size() == ps.n_neighbors(0, i));
        }

        std::cout << curParticlesPerDimension * curParticlesPerDimension * curParticlesPerDimension << ", " << durationBruteForce << ", " << durationCompactNSearch << std::endl;
    }
}

void investigateSmoothingLengths() {
    std::vector<std::array<double, 3>> particles = generateParticles(10, 1.0);
    const std::size_t particlesSize = particles.size();
    std::vector<double> analyticalNorms;
    for (auto& particle : particles) {
        analyticalNorms.push_back(Eigen::Vector3d(particle.data()).norm());
    }
    std::vector<double> interpolatedNorms(10*10*10);
    for (double radius = 0.1; radius < 2.0; radius += 0.1) {
        CompactNSearch::NeighborhoodSearch nsearch(radius);
        unsigned int pointSetId = nsearch.add_point_set(particles.front().data(), particlesSize);
        nsearch.z_sort();
        auto const& d = nsearch.point_set(pointSetId);
        d.sort_field(particles.data());
        d.sort_field(analyticalNorms.data());
        nsearch.find_neighbors();

        CompactNSearch::PointSet const& ps = nsearch.point_set(pointSetId);

        auto t1 = std::chrono::high_resolution_clock::now();
        for (std::size_t i = 0; i < particlesSize; i++) {
            double sumOfWeights = 0;
            double kernelSum = 0;
            std::array<double, 3> particle = particles[i];
            for (std::size_t j = 0; j < ps.n_neighbors(0, i); j++) {
                const unsigned int pid = ps.neighbor(0, i, j);
                std::array<double, 3> neighbor = particles[pid];
                // radius = 2h !!!
                double kernel = learnSPH::Kernel::kernel(Eigen::Vector3d(particle.data()), Eigen::Vector3d(neighbor.data()), 0.5 * radius);
                kernelSum += kernel;
                sumOfWeights += kernel * Eigen::Vector3d(neighbor.data()).norm();
            }
            interpolatedNorms[i] = sumOfWeights / kernelSum;
        }
        auto t2 = std::chrono::high_resolution_clock::now();

        auto durationInterpolation = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
        double minError = DBL_MAX;
        double maxError = 0;
        double totalError = 0;
        for (std::size_t idx = 0; idx < particlesSize; idx++) {
            double diff = abs(analyticalNorms[idx] - interpolatedNorms[idx]);
            totalError += diff;
            if (maxError < diff) {
                maxError = diff;
            }
            if (minError > diff) {
                minError = diff;
            }
        }
        double avgError = totalError / particlesSize;
        std::cout << radius << ", " << minError << ", " << maxError << ", " << avgError << ", " << durationInterpolation << std::endl;
    }
}

int main()
{
    investigateSmoothingLengths();

    std::cout << "Comparing Neighborhood Searches:" << std::endl;
    compareNeighborhoodSearches(15, 1.0, 0.1);

    std::cout << "Welcome to the learnSPH framework!!" << std::endl;
    std::cout << "Generating a sample scene..." << std::endl;
    // Generate particles
    std::vector<Eigen::Vector3d> particles;
    for (int i = 0; i < 20; i++) {
        for (int j = 0; j < 20; j++) {
            for (int k = 0; k < 20; k++) {
                double x = static_cast<double> (rand()) / static_cast <float> (RAND_MAX);
                double y = static_cast<double> (rand()) / static_cast <float> (RAND_MAX);
                double z = static_cast<double> (rand()) / static_cast <float> (RAND_MAX);
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
            particles_vector_data[particle_i] = (Eigen::Vector3d(0.5, 0.5, 1.0) - particles[particle_i]).normalized();
        }

        // Save output
        const std::string filename = "../res/example_" + std::to_string(time_step) + ".vtk";
        learnSPH::save_particles_to_vtk(filename, particles, particles_scalar_data, particles_vector_data);
    }

    std::cout << "completed!" << std::endl;
    std::cout << "The scene files have been saved in the folder `<build_folder>/res`. You can visualize them with Paraview." << std::endl;

    return 0;
}
