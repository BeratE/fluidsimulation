#include <stdlib.h>     // rand
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>    // std::max
#include <Eigen/Dense>
#include <omp.h>

int main()
{    
    const int NUM_STEPS = 100000;
    const int REQ_THREADS = 8;
    omp_set_num_threads(REQ_THREADS);

    double step_size = 1.0/(double)NUM_STEPS;
    double sum = 0.0;

    int nthreads;

    // Version one
    double tbegin = omp_get_wtime();

    #pragma omp parallel
    {
        const int NUM_THREADS = omp_get_num_threads();
        const int ID = omp_get_thread_num();

        if (ID == 0) nthreads = NUM_THREADS;

        double s = 0.0;
        for (int i = ID; i < NUM_STEPS; i += NUM_THREADS) {
            const double x = (i + 0.5) * step_size;
            s += (4.0) / (1.0 + x * x);
        }

        #pragma omp critical
        sum += s;
    }

    double pi = step_size * sum;
    
    double tend = omp_get_wtime();

    printf("Num Threads: %d\n", nthreads);
    printf("Pi = %f\n", pi);
    printf("DeltaT : %f\n", tend-tbegin);

    // Version two
    tbegin = omp_get_wtime();

    sum = 0.0;
    #pragma omp parallel for reduction(+: sum)
    for (int i = 0; i < NUM_STEPS; i++) {
        const double x = (i + 0.5) * step_size;
        sum += (4.0) / (1.0 + x * x);            
    }

    pi = step_size * sum;

    tend = omp_get_wtime();
    
    printf("Pi = %f\n", pi);
    printf("DeltaT : %f\n", tend-tbegin);

    return 0;
}
