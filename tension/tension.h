#pragma once
#include <Eigen/Dense>


namespace tension {
	// TODO: Find out appropriate values
	namespace Cohesion {
		namespace Kernel {
			double weight(const double r, const double c);
		}

		Eigen::Vector3d forceCohesion(const double mass_i, 
									  const double mass_j,
									  const Eigen::Vector3d pos_i,
									  const Eigen::Vector3d pos_j,
									  const double gamma,
									  const double c);
	}

	namespace Curvature {
		Eigen::Vector3d forceCurvature(const double mass_i,
			const Eigen::Vector3d normal_i,
			const Eigen::Vector3d normal_j,
			const double gamma);
	}

	namespace Adhesion {
		namespace Kernel {
			double weight(const double r, const double c);
		}

		Eigen::Vector3d forceAdhesion(const double mass_i,
			const double repVolume_k,
			const Eigen::Vector3d pos_i,
			const Eigen::Vector3d pos_k,
			const double beta,
			const double c);
	}

	Eigen::Vector3d forceTension(const double restDensity, 
		const double density_i,
		const double density_j,
		const Eigen::Vector3d cohesionForce,
		const Eigen::Vector3d curvatureForce);
}
