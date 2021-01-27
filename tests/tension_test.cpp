#include "catch.hpp"
#include "config.h"
#include "vtk_writer.h"
#include "learnSPH/solver_sph.h"
#include "learnSPH/solver_pbf.h"
#include "learnSPH/system/fluidsystem.h"
#include "learnSPH/system/emitter.h"

using namespace learnSPH;
using namespace learnSPH::System;

TEST_CASE("SPH Solver with Surface Tension", "[sph_tension]") {
	SECTION("SPH_TENSION") {
		std::cout << "Testing SPH with tension and adhesion" << std::endl;

		const double particleDiameter = 0.1;

		FluidSystem particles = Emitter().sampleFluidBox(Eigen::Vector3d(0, 0, 0),
														 Eigen::Vector3d(1, 1, 1),
														 particleDiameter);
		SolverSPH solver(particles);
		solver.setSnapShotAfterMS(40);
		solver.setParamStiffness(3000);
		solver.setFluidViscosity(0.1);
		solver.enableAdhesion(false);
		solver.enableGravity(false);
		solver.enableSmoothing(true);
		solver.enableTension(false);


		/*SECTION("SPH_NO_TENSION") {
			solver.newRun("new_sph_no_tension", 6000);
		}*/

		SECTION("SPH_TENSION") {
			solver.enableTension(true);
			solver.newRun("new_sph_tension", 1000);
		}
	}
}
//
//TEST_CASE("PBF Solver with Surface Tension", "[pbf_tension]") {
//	SECTION("PBF_TENSION") {
//		std::cout << "Testing PBF with tension and adhesion" << std::endl;
//
//		const double particleDiameter = 0.1;
//
//		FluidSystem particles = Emitter().sampleFluidBox(Eigen::Vector3d(0, 0, 0),
//			Eigen::Vector3d(1, 1, 1),
//			particleDiameter);
//		
//		SolverPBF solver(particles);
//
//		solver.setParamSmoothing(0.1);
//		solver.setSnapShotAfterMS(40);
//		solver.setFluidViscosity(0.0);
//		solver.enableAdhesion(false);
//		solver.enableGravity(false);
//		solver.enableSmoothing(true);
//		solver.enableTension(false);
//		solver.setNumIterations(3);
//
//
//		SECTION("PBF_NO_TENSION") {
//			solver.run("pbf_no_tension", 6000);
//		}
//
//		SECTION("PBF_NO_TENSION_NO_SMOOTHING") {
//			solver.enableSmoothing(false);
//			solver.run("pbf_no_tension_no_smoothing", 6000);
//		}
//
//		SECTION("PBF_TENSION") {
//			solver.enableTension(true);
//			solver.run("pbf_tension", 6000);
//		}
//
//		SECTION("PBF_TENSION_NO_SMOOTHING") {
//			solver.enableTension(true);
//			solver.enableSmoothing(false);
//			solver.run("pbf_tension_no_smoothing", 6000);
//		}
//	}
//}

//TEST_CASE("Tension on funnel", "[tension_funnel]") {
//	std::cout << "Funnel with tension and adhesion" << std::endl;
//
//	const double particleDiameter = 0.05;
//	const double particelDiameterBoundary = particleDiameter / 0.5;
//
//	FluidSystem particles = Emitter().sampleFluidBox(Eigen::Vector3d(0, 0, 0),
//		Eigen::Vector3d(1, 1, 1),
//		particleDiameter);
//	particles.setC(0.15); // Because particleDiameter went from 0.1 to 0.05
//
//	BoundarySystem plane1 = Emitter().sampleBoundaryPlane(Eigen::Vector3d(-1, 0, -1),
//		Eigen::Vector3d(0.3, -1, 0.3),
//		Eigen::Vector3d(2, 0, -1),
//		Eigen::Vector3d(0.7, -1, 0.3),
//		particelDiameterBoundary);
//	plane1.setViscosity(0.02);
//
//	BoundarySystem plane2 = Emitter().sampleBoundaryPlane(Eigen::Vector3d(2, 0, -1),
//		Eigen::Vector3d(0.7, -1, 0.3),
//		Eigen::Vector3d(2, 0, 2),
//		Eigen::Vector3d(0.7, -1, 0.7),
//		particelDiameterBoundary);
//	plane2.setViscosity(0.02);
//
//	BoundarySystem plane3 = Emitter().sampleBoundaryPlane(Eigen::Vector3d(2, 0, 2),
//		Eigen::Vector3d(0.7, -1, 0.7),
//		Eigen::Vector3d(-1, 0, 2),
//		Eigen::Vector3d(0.3, -1, 0.7),
//		particelDiameterBoundary);
//	plane3.setViscosity(0.02);
//
//	BoundarySystem plane4 = Emitter().sampleBoundaryPlane(Eigen::Vector3d(-1, 0, 2),
//		Eigen::Vector3d(0.3, -1, 0.7),
//		Eigen::Vector3d(-1, 0, -1),
//		Eigen::Vector3d(0.3, -1, 0.3),
//		particelDiameterBoundary);
//	plane4.setViscosity(0.02);
//
//	SECTION("SPH") {
//		SolverSPH solver(particles);
//		solver.addBoundary(plane1);
//		solver.addBoundary(plane2);
//		solver.addBoundary(plane3);
//		solver.addBoundary(plane4);
//
//		solver.setSnapShotAfterMS(40);
//		solver.setParamStiffness(3000);
//		solver.setFluidViscosity(0.1);
//		solver.enableAdhesion(false);
//		solver.enableGravity(true);
//		solver.enableSmoothing(true);
//		solver.enableTension(false);
//
//		SECTION("OFF") {
//			solver.run("sph_funnel_no_tension_no_adhesion", 4000);
//		}
//
//		SECTION("BETA250") {
//			solver.enableAdhesion(true);
//			solver.enableTension(true);
//			plane1.setBeta(250);
//			plane2.setBeta(250);
//			plane3.setBeta(250);
//			plane4.setBeta(250);
//			solver.run("sph_funnel_250_beta", 4000);
//		}
//
//		SECTION("BETA500") {
//			solver.enableAdhesion(true);
//			solver.enableTension(true);
//			plane1.setBeta(500);
//			plane2.setBeta(500);
//			plane3.setBeta(500);
//			plane4.setBeta(500);
//			solver.run("sph_funnel_500_beta", 4000);
//		}
//	}
//
//	SECTION("PBF") {
//		SolverPBF solver(particles);
//		solver.addBoundary(plane1);
//		solver.addBoundary(plane2);
//		solver.addBoundary(plane3);
//		solver.addBoundary(plane4);
//
//		solver.setParamSmoothing(0.25);
//		solver.setSnapShotAfterMS(40);
//		solver.setFluidViscosity(0.0);
//		solver.enableAdhesion(false);
//		solver.enableGravity(true);
//		solver.enableSmoothing(true);
//		solver.enableTension(false);
//		solver.setNumIterations(3);
//
//		SECTION("OFF") {
//			solver.run("pbf_funnel_no_tension_no_adhesion", 4000);
//		}
//
//		SECTION("BETA250") {
//			solver.enableAdhesion(true);
//			solver.enableTension(true);
//			plane1.setBeta(250);
//			plane2.setBeta(250);
//			plane3.setBeta(250);
//			plane4.setBeta(250);
//			solver.run("pbf_funnel", 4000);
//		}
//
//		SECTION("BETA500") {
//			solver.enableAdhesion(true);
//			solver.enableTension(true);
//			plane1.setBeta(500);
//			plane2.setBeta(500);
//			plane3.setBeta(500);
//			plane4.setBeta(500);
//			solver.run("pbf_funnel", 4000);
//		}
//	}
//
//}

//TEST_CASE("Surface Tension Complex Scene", "[tension_complex]") {
//	std::cout << "Complex scene with tension and adhesion" << std::endl;
//
//	const double particleDiameter = 0.05;
//	const double particelDiameterBoundary = particleDiameter / 0.5;
//
//	FluidSystem particles = Emitter().sampleFluidBox( Eigen::Vector3d(0.1, 1.1, 1.9),//Eigen::Vector3d(0.05, 1.05, 1.95),
//		Eigen::Vector3d(1.9, 1.9, 0.1),
//		particleDiameter);
//	particles.setC(0.15);
//
//	BoundarySystem bigBoundary = Emitter().sampleBoundaryHollowBox(Eigen::Vector3d(0.0, 0.0, 0.0),
//		Eigen::Vector3d(4.0, 2.0, 2.0),
//		particelDiameterBoundary);
//	bigBoundary.setViscosity(0.02);
//
//	BoundarySystem plane1 = Emitter().sampleBoundaryPlane(Eigen::Vector3d(0.0, 1.0, 2.0),
//		Eigen::Vector3d(2.0, 1.0, 2.0),
//		Eigen::Vector3d(0.0, 1.0, 0.0),
//		Eigen::Vector3d(2.0, 1.0, 0.0),
//		particelDiameterBoundary);
//	plane1.setViscosity(0.02);
//
//	BoundarySystem plane2 = Emitter().sampleBoundaryPlane(Eigen::Vector3d(2.0, 2.0, 2.0),
//		Eigen::Vector3d(2.0, 1.5, 2.0),
//		Eigen::Vector3d(2.0, 2.0, 0.0),
//		Eigen::Vector3d(2.0, 1.5, 0.0),
//		particelDiameterBoundary);
//	plane2.setViscosity(0.02);
//
//	BoundarySystem plane3 = Emitter().sampleBoundaryPlane(Eigen::Vector3d(2.0, 1.5, 2.0),
//		Eigen::Vector3d(2.0, 1.0, 2.0),
//		Eigen::Vector3d(2.0, 1.5, 1.25),
//		Eigen::Vector3d(2.0, 1.0, 1.25),
//		particelDiameterBoundary);
//	plane3.setViscosity(0.02);
//
//	BoundarySystem plane4 = Emitter().sampleBoundaryPlane(Eigen::Vector3d(2.0, 1.5, 0.75),
//		Eigen::Vector3d(2.0, 1.0, 0.75),
//		Eigen::Vector3d(2.0, 1.5, 0.0),
//		Eigen::Vector3d(2.0, 1.0, 0.0),
//		particelDiameterBoundary);
//	plane4.setViscosity(0.02);
//
//	std::stringstream filepath;
//	filepath << SOURCE_DIR << "/res/" << "armadillo_low.obj";
//	BoundarySystem armadillo = Emitter().sampleBoundaryMesh(filepath.str(), particelDiameterBoundary);
//	armadillo.setViscosity(0.08);
//
//	Eigen::Matrix4d transform = Eigen::Matrix4d::Identity();
//	transform(0, 3) = 2.5;
//	transform(1, 3) = 0.1;
//	transform(2, 3) = 1.0;
//	armadillo.transform(transform);
//
//	SECTION("SPH") {
//		SolverSPH solver(particles);
//		solver.addBoundary(bigBoundary);
//		solver.addBoundary(plane1);
//		solver.addBoundary(plane2);
//		solver.addBoundary(plane3);
//		solver.addBoundary(plane4);
//		solver.addBoundary(armadillo);
//
//		solver.setSnapShotAfterMS(40);
//		solver.setParamStiffness(3000);
//		solver.setFluidViscosity(0.1);
//		solver.enableAdhesion(false);
//		solver.enableGravity(true);
//		solver.enableSmoothing(true);
//		solver.enableTension(false);
//
//		SECTION("SPH_tension_complex_no_tension") {
//			solver.run("sph_tension_complex_no_tension", 6000);
//		}
//
//		SECTION("SPH_tension_complex_no_tension") {
//			solver.enableAdhesion(true);
//			solver.enableTension(true);
//			bigBoundary.setBeta(500);
//			plane1.setBeta(500);
//			plane2.setBeta(500);
//			plane3.setBeta(500);
//			plane4.setBeta(500);
//			solver.run("sph_tension_complex", 6000);
//		}
//	}
//
//
//	SECTION("PBF") {
//		SolverPBF solver(particles);
//		solver.addBoundary(bigBoundary);
//		solver.addBoundary(plane1);
//		solver.addBoundary(plane2);
//		solver.addBoundary(plane3);
//		solver.addBoundary(plane4);
//		solver.addBoundary(armadillo);
//
//		solver.setParamSmoothing(0.25);
//		solver.setSnapShotAfterMS(40);
//		solver.setFluidViscosity(0.0);
//		solver.enableAdhesion(false);
//		solver.enableGravity(true);
//		solver.enableSmoothing(true);
//		solver.enableTension(false);
//		solver.setNumIterations(3);
//
//		SECTION("PBF_tension_complex") {
//			solver.run("pbf_tension_complex_no_tension", 6000);
//		}
//
//		SECTION("SPH_tension_complex_no_tension") {
//			solver.enableAdhesion(true);
//			solver.enableTension(true);
//			bigBoundary.setBeta(500);
//			plane1.setBeta(500);
//			plane2.setBeta(500);
//			plane3.setBeta(500);
//			plane4.setBeta(500);
//			solver.run("pbf_tension_complex", 6000);
//		}
//	}
//}
