#include <Eigen/Dense>
#include "motionIntegrator.h"
#include "kernel.h"

using namespace learnSPH;
using namespace learnSPH::MotionIntegrator;

void
learnSPH::MotionIntegrator::semiImplicitEulerIntegrator(ParticleSystem::FluidSystem& fluid, const double minimumTimeStep, const CompactNSearch::NeighborhoodSearch& nsearch, double epsilon) {
	// smoothin length 
	const double h = Kernel::Parameter::TUNING * fluid.particleRadius * 2;
	const double timeStepCFL = fluid.getTimeCFL();
	const double timeStep = minimumTimeStep < timeStepCFL ? minimumTimeStep : timeStepCFL;

	// get neighborhood information of fluid particle point set
	CompactNSearch::PointSet const& fluidPS = nsearch.point_set(fluid.id);
	fluid.positions.resize(fluidPS.n_points());
	fluid.velocities.resize(fluidPS.n_points());
	fluid.accelerations.resize(fluidPS.n_points());

	for (size_t fpI = 0; fpI < fluidPS.n_points(); fpI++) {
		fluid.velocities[fpI] = fluid.velocities[fpI] + timeStep * fluid.accelerations[fpI];
	}

	// iterate fluid particles 
	for (size_t fpI = 0; fpI < fluidPS.n_points(); fpI++) {
		Eigen::Vector3d fpPos = fluid.positions[fpI];
		Eigen::Vector3d fpVelStar = Eigen::Vector3d::Zero();
		fpVelStar += fluid.velocities[fpI];
		Eigen::Vector3d fpVelSumOverNeighbors = Eigen::Vector3d::Zero();
		/* 
		 * No need to pay attention to the particle itself when calculating fpSmoothedVelSum,
		 * because vj - vi == 0 if i == j
		 */
		for (size_t fpN = 0; fpN < fluidPS.n_neighbors(fluid.id, fpI); fpN++) {
			const unsigned int fnI = fluidPS.neighbor(fluid.id, fpI, fpN);
			fpVelSumOverNeighbors += 2 * fluid.particleMass * pow(fluid.densities[fpI] + fluid.densities[fnI], -1) * Kernel::CubicSpline::weight(fpPos, fluid.positions[fnI], h) * (fluid.velocities[fnI] - fluid.velocities[fpI]);
		}

		fpVelSumOverNeighbors *= epsilon;
		fpVelStar += fpVelSumOverNeighbors;

		const Eigen::Vector3d newFpPos = fpPos + timeStep * fpVelStar;
		fluid.positions[fpI] = newFpPos;
		fluid.velocities[fpI] = fluid.velocities[fpI];
	}
}