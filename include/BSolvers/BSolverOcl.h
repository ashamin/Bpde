#ifndef BSOLVER_OCL_H
#define BSOLVER_OCL_H

#include "BArea.h"
#include "BSolver.h"

#include "CL/cl.h"
#include "omp.h"

#include <iostream>
#include "math.h"

namespace Bpde {

class BSolverOcl : public BSolver
{
public:
    BSolverOcl(const BArea& area, cl_device_id deviceId, int threadsNum);
    virtual ~BSolverOcl();

    double* solve();

    virtual double exec_time();
    int it_num();
private:
    inline void prepareIteration();
    int iterations;
    double time;
    int threadsNum;
    BArea area;

    double *H, *Ha;
    double *x, *y;
    double *V;
    double *mu;
    double *dx_d, *dx_l, *dx_u, *dy_d, *dy_l, *dy_u;
    double *loc_c, *loc_d, *b;
    int I, J, T;
    double dt;
    int n;
    double t;

    cl_device_id deviceId;
    cl_context context;
    cl_command_queue commandQueue;
    cl_program program;
    cl_kernel kernel;
};

} // namespace Bpde

#endif
