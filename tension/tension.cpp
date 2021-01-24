#define _USE_MATH_DEFINES
#include "tension.h"
#include <math.h>

using namespace tension;

double Cohesion::Kernel::weight(const double r) {
	const double alpha = 32.0 / (M_PI * pow(Parameter::C, 9));
	if (0.0 <= r && r <= Parameter::C / 2.0) {
		return alpha * 2 * pow((Parameter::C - r), 3) * pow(r, 3) - pow(Parameter::C, 6) / 64.0;
	}
	else if (Parameter::C / 2.0 < r && r <= Parameter::C) {
		return alpha * pow((Parameter::C - r), 3) * pow(r, 3);
	}
	return 0.0;
}

Eigen::Vector3d Cohesion::forceCohesion(const double mass_i,
	const double mass_j,
	const Eigen::Vector3d pos_i,
	const Eigen::Vector3d pos_j) {

	const Eigen::Vector3d diff = pos_i - pos_j;
	const double diffNorm = diff.norm();

	return -Parameter::GAMMA * mass_i * mass_i * Cohesion::Kernel::weight(diffNorm) * diff.normalized();
}

Eigen::Vector3d Curvature::forceCurvature(const double mass_i,
	const Eigen::Vector3d normal_i,
	const Eigen::Vector3d normal_j) {
	return -Parameter::GAMMA * mass_i * (normal_i - normal_j);
}

double Adhesion::Kernel::weight(const double r) {
	const double alpha = 0.007 / pow(Parameter::C, 3.25);

	if (Parameter::C / 2.0 <= r && r <= Parameter::C) {
		return alpha * pow((-4.0 * pow(r, 2)) / Parameter::C + 6.0 * r - 2.0 * Parameter::C, 1.0 / 4.0);
	}
	return 0;
}

Eigen::Vector3d Adhesion::forceAdhesion(const double mass_i,
	const double repVolume_k,
	const Eigen::Vector3d pos_i,
	const Eigen::Vector3d pos_k) {

	const Eigen::Vector3d diff = pos_i - pos_k;
	const double diffNorm = diff.norm();

	return -Parameter::BETA * mass_i * repVolume_k * Kernel::weight(diffNorm) * diff.normalized();
}

Eigen::Vector3d tension::forceTension(const double restDensity,
	const double density_i,
	const double density_j,
	const Eigen::Vector3d cohesionForce,
	const Eigen::Vector3d curvatureForce) {
	
	const double K = (2.0 * restDensity) / (density_i + density_j);
	return K * (cohesionForce + curvatureForce);
}