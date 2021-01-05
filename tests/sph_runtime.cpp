#include "catch.hpp"
#include <omp.h>
#include <iostream>
#include "config.h"
#include "vtk_writer.h"
#include "learnSPH/kernel.h"
#include "learnSPH/solver_sph.h"
#include "learnSPH/system/emitter.h"
#include "learnSPH/system/particlesystem.h"

TEST_CASE("RUN_TIME", "[runtime]")
{
    const double start_t = omp_get_wtime();


    const double end_t = omp_get_wtime();
    const double delta_t = end_t - start_t;

    std::cout << "Runtime: " << delta_t << std::endl;
}
