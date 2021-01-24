#pragma once
#include <Eigen/Dense>


namespace tension {
	// TODO: Find out appropriate values
	namespace Parameter {
		constexpr double GAMMA = 0.0;
		constexpr double C = 0.0;
		constexpr double BETA = 0.0;
	}
	namespace Cohesion {
		namespace Kernel {
			double weight(const double r);
		}

		Eigen::Vector3d forceCohesion(const double mass_i, 
									  const double mass_j,
									  const Eigen::Vector3d pos_i,
									  const Eigen::Vector3d pos_j);
	}

	namespace Curvature {
		Eigen::Vector3d forceCurvature(const double mass_i,
			const Eigen::Vector3d normal_i,
			const Eigen::Vector3d normal_j);
	}

	namespace Adhesion {
		namespace Kernel {
			double weight(const double r);
		}

		Eigen::Vector3d forceAdhesion(const double mass_i,
			const double repVolume_k,
			const Eigen::Vector3d pos_i,
			const Eigen::Vector3d pos_k);
	}

	Eigen::Vector3d forceTension(const double restDensity, 
		const double density_i,
		const double density_j,
		const Eigen::Vector3d cohesionForce,
		const Eigen::Vector3d curvatureForce);
}
