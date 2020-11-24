#pragma once
#include "particlesystem.h"
#include <CompactNSearch/CompactNSearch.h>


namespace learnSPH {
	namespace Simulator {
		/*
		 * Estimate the new motion state of all the particles of a fluid using the semi-implicit Euler.
		 * @param &fluid - Fluid for whose particles the motion update should be calculated.
		 * @param minimumTimeStep - Lower bound for the time step that should be used in the integration.
		 * @param &nsearch - Compact neighborhood information.
		 * @param epsilon - Scaling factor for velocity estimation in XSPH. 
		 */
		void semiImplicitEuler(ParticleSystem::FluidSystem& fluid, const double minimumTimeStep, const CompactNSearch::NeighborhoodSearch& nsearch, double epsilon);
	}
}