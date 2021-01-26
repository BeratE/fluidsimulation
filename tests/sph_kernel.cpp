
#include "catch.hpp"
#include <tuple>
#include <stdio.h>
#include <Eigen/Dense>
#include "learnSPH/kernel.h"

using namespace learnSPH::Kernel;

//TEST_CASE("Cubic", "[kernel]")
//{
//    SECTION("CubicKernel") {
//        const double h = 0.3;
//        SECTION("ApproxInt", "[Approximate integral of s(q) using Simpson rule]") {
//            // Approximate integral using simpson rule
//            std::vector<std::tuple<double, double>> values;
//            for (double dist = -3.0; dist < 3.0; dist += 0.1) {
//                values.push_back(std::make_tuple(dist, CubicSpline::cubicSpline(fabs(dist))));
//            }
//
//            double integral = 0.0;
//            for (int i = 0; i < values.size()-1; i++) {
//                auto [dist_a, val_a] = values[i];
//                auto [dist_b, val_b] = values[i+1];
//                double m = (dist_b + dist_a)/2.0;
//                double fm = CubicSpline::cubicSpline(fabs(m));
//                integral += ((dist_b - dist_a)/6) * (val_a + 4*fm + val_b);
//            }
//
//            //REQUIRE (integral == Approx(1.0).epsilon(0.05));
//        }
//    
//        SECTION("Manual", "[Compare to manually calculated values]") {
//            REQUIRE(CubicSpline::weight
//                    (Eigen::Vector3d(0.0, 0.0, 0.0),
//                     Eigen::Vector3d(2.5*h, 0.0, 0.0), // >2h
//                     h) == 0);
//            REQUIRE(CubicSpline::weight
//                    (Eigen::Vector3d(0.0, 0.0, 0.0),
//                     Eigen::Vector3d(0.5*h, 0.0, 0.0), // <1h
//                     h) == Approx(0.2287852307 * pow(h,-3)));
//            REQUIRE(CubicSpline::weight
//                    (Eigen::Vector3d(0.0, 0.0, 0.0),
//                     Eigen::Vector3d(1.5*h, 0.0, 0.0), // 1 < h < 2
//                     h) == Approx(0.009947183943* pow(h,-3)));
//        }
//    
//        SECTION("Properties", "[Kernel properties at random points]") {
//            const int NUM_PARTICLES = 100;
//            // Generate particles in a cube of size 1x1x1
//            std::vector<Eigen::Vector3d> particles;
//            for (int i = 0; i < NUM_PARTICLES; i++) {
//                particles.push_back(Eigen::Vector3d(rand()/(double)RAND_MAX,
//                                                    rand()/(double)RAND_MAX,
//                                                    rand()/(double)RAND_MAX));
//      
//            }
//        
//            // Compact support of cubic function is 2h
//            double alpha = 2*h;
//            double sum = 0.0;
//            for (int i = 0; i < NUM_PARTICLES-1; i++) {
//                for (int j = i + 1; j < NUM_PARTICLES; j++) {
//                    double k1 = CubicSpline::weight(particles[i], particles[j], h);
//                    double k2 = CubicSpline::weight(particles[j], particles[i], h);
//                    // Symmetry
//                    REQUIRE(k1 == k2);
//                    // Non-negative
//                    REQUIRE(k1 >= 0);
//                    // Compact support
//                    if ((particles[i] - particles[j]).norm() > alpha)
//                        REQUIRE(k1 == 0);
//
//                    sum += k1;
//                }
//            }
//        }
//    }
//}
//
//Eigen::Vector3d cubicKernelGradApprox(Eigen::Vector3d x_i, Eigen::Vector3d x_j,
//                                      double h, double eps)
//{
//    Eigen::Vector3d e_x(eps, 0.0, 0.0);
//    Eigen::Vector3d e_y(0.0, eps, 0.0);
//    Eigen::Vector3d e_z(0.0, 0.0, eps);
//    Eigen::Vector3d grad;
//    grad[0] = CubicSpline::weight(x_i + e_x, x_j, h) - CubicSpline::weight(x_i - e_x, x_j, h);
//    grad[1] = CubicSpline::weight(x_i + e_y, x_j, h) - CubicSpline::weight(x_i - e_y, x_j, h);
//    grad[2] = CubicSpline::weight(x_i + e_z, x_j, h) - CubicSpline::weight(x_i - e_z, x_j, h);
//    grad /= (2*eps);
//    return grad;
//}
//
//TEST_CASE("CubicGradient", "[kernel][gradient]") {
//    SECTION("CubicKernelGradient") {
//        const double eps = pow(10, -6);
//        const double h = 0.4;
//
//        const int NUM_PARTICLES = 20;
//        // Generate particles in a cube of size 1x1x1
//        std::vector<Eigen::Vector3d> particles;
//        for (int i = 0; i < NUM_PARTICLES; i++) {
//            particles.push_back(Eigen::Vector3d(rand()/(double)RAND_MAX,
//                                                rand()/(double)RAND_MAX,
//                                                rand()/(double)RAND_MAX));
//        }
//
//        SECTION("CompGradApprox", "Compare analytical solution with approximation at random points.") {
//            for (int i = 0; i < NUM_PARTICLES-1; i++) {
//                for (int j = i + 1; j < NUM_PARTICLES; j++) {
//                    Eigen::Vector3d gradApprox = cubicKernelGradApprox(particles[i],
//                                                                       particles[j], h, eps);
//                    Eigen::Vector3d gradCubic = CubicSpline::gradWeight(particles[i],
//                                                                        particles[j], h);
//
//                    for (int k = 0; k < 3; k++) {
//                        INFO("q = " << (particles[i] - particles[j]).norm() / h
//                             << " Grad = " << gradCubic[k] << " Approx~ "
//                             << gradApprox[k] << " Diff: "
//                             << fabs(gradCubic[k] - gradApprox[k]) << "\n");
//
//                        REQUIRE(gradCubic[k] == Approx(gradApprox[k]).epsilon(0.01));
//                    }
//                }
//            }
//        }
//    }
//}
