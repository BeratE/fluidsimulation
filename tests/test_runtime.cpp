#include "catch.hpp"
#include "kernel.h"
#include "solver.h"
#include "system/emitter.h"
#include "system/particlesystem.h"
#include "util/config.h"
#include "util/vtk_writer.h"
#include <omp.h>
#include <iostream>

TEST_CASE("RUN_TIME", "[runtime]")
{
    const double start_t = omp_get_wtime();


    const double end_t = omp_get_wtime();
    const double delta_t = end_t - start_t;

    std::cout << "Runtime: " << delta_t << std::endl;
}
