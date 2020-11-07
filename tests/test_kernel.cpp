#include "catch.hpp"
#include <Eigen/Dense>
#include "kernel.h"

TEST_CASE( "Cubic kernel function", "[kernel]" )
{

    const double h = 0.05;
    const int NUM_PARTICLES = 30;
    // Generate particles in a cube of size 1x1x1
    std::vector<Eigen::Vector3d> particles;
    for (int i = 0; i < NUM_PARTICLES; i++) {
        particles.push_back(Eigen::Vector3d(rand()/(double)RAND_MAX,
                                            rand()/(double)RAND_MAX,
                                            rand()/(double)RAND_MAX));
      
    }
    
    SECTION("Compact support of kernel") {
        // Compact support of cubic function is 2h
        double alpha = 2*h;
        for (int i = 0; i < NUM_PARTICLES; i++) {
            for (int j = i+1; j < NUM_PARTICLES; j++) {
                if ((particles[i] - particles[j]).norm() > alpha) 
                    REQUIRE(learnSPH::Kernel::kernel(particles[i], particles[j], h) == 0);
          }
        }
    }
}
