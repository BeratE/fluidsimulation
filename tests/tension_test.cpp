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


		SECTION("SPH_NO_TENSION") {
			solver.newRun("sph_no_tension", 1000);
		}

		SECTION("SPH_TENSION") {
			solver.enableTension(true);
			solver.newRun("sph_tension", 1000);
		}
	}
}

TEST_CASE("PBF Solver with Surface Tension", "[pbf_tension]") {
	SECTION("PBF_TENSION") {
		std::cout << "Testing PBF with tension and adhesion" << std::endl;

		const double particleDiameter = 0.1;

		FluidSystem particles = Emitter().sampleFluidBox(Eigen::Vector3d(0, 0, 0),
			Eigen::Vector3d(1, 1, 1),
			particleDiameter);
		
		SolverPBF solver(particles);

		solver.setParamSmoothing(0.1);
		solver.setSnapShotAfterMS(40);
		solver.setFluidViscosity(0.0);
		solver.enableAdhesion(false);
		solver.enableGravity(false);
		solver.enableSmoothing(true);
		solver.enableTension(false);
		solver.setNumIterations(3);


		SECTION("PBF_NO_TENSION") {
			solver.newRun("pbf_no_tension", 1000);
		}

		SECTION("PBF_TENSION") {
			solver.enableTension(true);
			solver.newRun("pbf_tension", 1000);
		}
	}
}

TEST_CASE("Tension on funnel", "[tension_funnel]") {
	std::cout << "Funnel with tension and adhesion" << std::endl;

	const double particleDiameter = 0.05;
	const double particelDiameterBoundary = particleDiameter / 0.5;

	FluidSystem particles = Emitter().sampleFluidBox(Eigen::Vector3d(0, 0, 0),
		Eigen::Vector3d(1, 1, 1),
		particleDiameter);
	particles.setGamma(0.2);

	BoundarySystem plane1 = Emitter().sampleBoundaryPlane(Eigen::Vector3d(-1, 0, -1),
		Eigen::Vector3d(0.3, -1, 0.3),
		Eigen::Vector3d(2, 0, -1),
		Eigen::Vector3d(0.7, -1, 0.3),
		particelDiameterBoundary);
	plane1.setViscosity(0.02);

	BoundarySystem plane2 = Emitter().sampleBoundaryPlane(Eigen::Vector3d(2, 0, -1),
		Eigen::Vector3d(0.7, -1, 0.3),
		Eigen::Vector3d(2, 0, 2),
		Eigen::Vector3d(0.7, -1, 0.7),
		particelDiameterBoundary);
	plane2.setViscosity(0.02);

	BoundarySystem plane3 = Emitter().sampleBoundaryPlane(Eigen::Vector3d(2, 0, 2),
		Eigen::Vector3d(0.7, -1, 0.7),
		Eigen::Vector3d(-1, 0, 2),
		Eigen::Vector3d(0.3, -1, 0.7),
		particelDiameterBoundary);
	plane3.setViscosity(0.02);

	BoundarySystem plane4 = Emitter().sampleBoundaryPlane(Eigen::Vector3d(-1, 0, 2),
		Eigen::Vector3d(0.3, -1, 0.7),
		Eigen::Vector3d(-1, 0, -1),
		Eigen::Vector3d(0.3, -1, 0.3),
		particelDiameterBoundary);
	plane4.setViscosity(0.02);

	BoundarySystem box = Emitter().sampleBoundaryHollowBox(Eigen::Vector3d(-1, -1.5, 2),
		Eigen::Vector3d(2, 1.2, -1),
		particelDiameterBoundary);
	box.setViscosity(0.02);

	SECTION("SPH") {
		SolverSPH solver(particles);
		solver.addBoundary(plane1);
		solver.addBoundary(plane2);
		solver.addBoundary(plane3);
		solver.addBoundary(plane4);
		solver.addBoundary(box);

		solver.setSnapShotAfterMS(40);
		solver.setParamStiffness(3000);
		solver.setFluidViscosity(0.1);
		solver.enableAdhesion(false);
		solver.enableGravity(true);
		solver.enableSmoothing(true);
		solver.enableTension(false);

		SECTION("OFF") {
			solver.newRun("sph_funnel_no_tension_no_adhesion", 3000);
		}

		SECTION("BETA250") {
			solver.enableAdhesion(true);
			solver.enableTension(true);
			solver.newRun("sph_funnel_1_beta", 3000);
		}
	}

	SECTION("PBF") {
		SolverPBF solver(particles);
		solver.addBoundary(plane1);
		solver.addBoundary(plane2);
		solver.addBoundary(plane3);
		solver.addBoundary(plane4);
		solver.addBoundary(box);

		solver.setParamSmoothing(0.25);
		solver.setSnapShotAfterMS(40);
		solver.setFluidViscosity(0.0);
		solver.enableAdhesion(false);
		solver.enableGravity(true);
		solver.enableSmoothing(true);
		solver.enableTension(false);
		solver.setNumIterations(3);

		SECTION("OFF") {
			solver.newRun("pbf_funnel_no_tension_no_adhesion", 3000);
		}

		SECTION("BETA250") {
			solver.enableAdhesion(true);
			solver.enableTension(true);
			solver.newRun("pbf_funnel_1_beta", 3000);
		}
	}

}

TEST_CASE("Surface Tension Complex Scene", "[tension_complex]") {
	std::cout << "Complex scene with tension and adhesion" << std::endl;

	const double particleDiameter = 0.05;
	const double particelDiameterBoundary = particleDiameter / 0.5;

	FluidSystem particles = Emitter().sampleFluidBox( Eigen::Vector3d(0.1, 1.1, 1.9),//Eigen::Vector3d(0.05, 1.05, 1.95),
		Eigen::Vector3d(1.9, 1.9, 0.1),
		particleDiameter);
	particles.setGamma(0.2);

	BoundarySystem bigBoundary = Emitter().sampleBoundaryHollowBox(Eigen::Vector3d(0.0, 0.0, 0.0),
		Eigen::Vector3d(4.0, 2.0, 2.0),
		particelDiameterBoundary);
	bigBoundary.setViscosity(0.02);

	BoundarySystem plane1 = Emitter().sampleBoundaryPlane(Eigen::Vector3d(0.0, 1.0, 2.0),
		Eigen::Vector3d(2.0, 1.0, 2.0),
		Eigen::Vector3d(0.0, 1.0, 0.0),
		Eigen::Vector3d(2.0, 1.0, 0.0),
		particelDiameterBoundary);
	plane1.setViscosity(0.02);

	BoundarySystem plane2 = Emitter().sampleBoundaryPlane(Eigen::Vector3d(2.0, 2.0, 2.0),
		Eigen::Vector3d(2.0, 1.5, 2.0),
		Eigen::Vector3d(2.0, 2.0, 0.0),
		Eigen::Vector3d(2.0, 1.5, 0.0),
		particelDiameterBoundary);
	plane2.setViscosity(0.02);

	BoundarySystem plane3 = Emitter().sampleBoundaryPlane(Eigen::Vector3d(2.0, 1.5, 2.0),
		Eigen::Vector3d(2.0, 1.0, 2.0),
		Eigen::Vector3d(2.0, 1.5, 1.25),
		Eigen::Vector3d(2.0, 1.0, 1.25),
		particelDiameterBoundary);
	plane3.setViscosity(0.02);

	BoundarySystem plane4 = Emitter().sampleBoundaryPlane(Eigen::Vector3d(2.0, 1.5, 0.75),
		Eigen::Vector3d(2.0, 1.0, 0.75),
		Eigen::Vector3d(2.0, 1.5, 0.0),
		Eigen::Vector3d(2.0, 1.0, 0.0),
		particelDiameterBoundary);
	plane4.setViscosity(0.02);

	std::stringstream filepath;
	filepath << SOURCE_DIR << "/res/" << "armadillo_low.obj";
	BoundarySystem armadillo = Emitter().sampleBoundaryMesh(filepath.str(), particelDiameterBoundary);
	armadillo.setViscosity(0.08);

	Eigen::Matrix4d transform = Eigen::Matrix4d::Identity();
	transform(0, 3) = 2.5;
	transform(1, 3) = 0.1;
	transform(2, 3) = 1.0;
	armadillo.transform(transform);

	SECTION("SPH") {
		SolverSPH solver(particles);
		solver.addBoundary(bigBoundary);
		solver.addBoundary(plane1);
		solver.addBoundary(plane2);
		solver.addBoundary(plane3);
		solver.addBoundary(plane4);
		solver.addBoundary(armadillo);

		solver.setSnapShotAfterMS(40);
		solver.setParamStiffness(3000);
		solver.setFluidViscosity(0.1);
		solver.enableAdhesion(false);
		solver.enableGravity(true);
		solver.enableSmoothing(true);
		solver.enableTension(false);

		SECTION("SPH_tension_complex_no_tension") {
			solver.newRun("sph_tension_complex_no_tension", 3000);
		}

		SECTION("SPH_tension_complex_no_tension") {
			solver.enableAdhesion(true);
			solver.enableTension(true);
			solver.newRun("sph_tension_complex", 3000);
		}
	}


	SECTION("PBF") {
		SolverPBF solver(particles);
		solver.addBoundary(bigBoundary);
		solver.addBoundary(plane1);
		solver.addBoundary(plane2);
		solver.addBoundary(plane3);
		solver.addBoundary(plane4);
		solver.addBoundary(armadillo);

		solver.setParamSmoothing(0.25);
		solver.setSnapShotAfterMS(40);
		solver.setFluidViscosity(0.0);
		solver.enableAdhesion(false);
		solver.enableGravity(true);
		solver.enableSmoothing(true);
		solver.enableTension(false);
		solver.setNumIterations(3);

		SECTION("PBF_tension_complex") {
			solver.newRun("pbf_tension_complex_no_tension", 3000);
		}

		SECTION("SPH_tension_complex_no_tension") {
			solver.enableAdhesion(true);
			solver.enableTension(true);
			solver.newRun("pbf_tension_complex", 3000);
		}
	}
}
