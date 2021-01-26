#define _USE_MATH_DEFINES
#include "tension.h"
#include <math.h>

using namespace tension;

double Cohesion::Kernel::weight(const double r, const double c) {
	const double alpha = 32.0 / (M_PI * pow(c, 9));
	if (0.0 <= r && r <= c / 2.0) {
		return alpha * 2 * pow((c - r), 3) * pow(r, 3) - pow(c, 6) / 64.0;
	}
	else if (c / 2.0 < r && r <= c) {
		return alpha * pow((c - r), 3) * pow(r, 3);
	}
	return 0.0;
}

Eigen::Vector3d Cohesion::forceCohesion(const double mass_i,
	const double mass_j,
	const Eigen::Vector3d pos_i,
	const Eigen::Vector3d pos_j,
	const double gamma,
	const double c) {

	const Eigen::Vector3d diff = pos_i - pos_j;
	const double diffNorm = diff.norm();

	return -gamma * mass_i * mass_i * Cohesion::Kernel::weight(diffNorm, c) * diff.normalized();
}

Eigen::Vector3d Curvature::forceCurvature(const double mass_i,
	const Eigen::Vector3d normal_i,
	const Eigen::Vector3d normal_j,
	const double gamma) {
	return -gamma * mass_i * (normal_i - normal_j);
}

double Adhesion::Kernel::weight(const double r, const double c) {
	const double alpha = 0.007 / pow(c, 3.25);

	if (c / 2.0 <= r && r <= c) {
		return alpha * pow((-4.0 * pow(r, 2)) / c + 6.0 * r - 2.0 * c, 1.0 / 4.0);
	}
	return 0;
}

Eigen::Vector3d Adhesion::forceAdhesion(const double mass_i,
	const double repVolume_k,
	const Eigen::Vector3d pos_i,
	const Eigen::Vector3d pos_k,
	const double beta,
	const double c) {

	const Eigen::Vector3d diff = pos_i - pos_k;
	const double diffNorm = diff.norm();

	return -beta * mass_i * repVolume_k * Kernel::weight(diffNorm, c) * diff.normalized();
}

Eigen::Vector3d tension::forceTension(const double restDensity,
	const double density_i,
	const double density_j,
	const Eigen::Vector3d cohesionForce,
	const Eigen::Vector3d curvatureForce) {
	
	const double K = (2.0 * restDensity) / (density_i + density_j);
	return K * (cohesionForce + curvatureForce);
}