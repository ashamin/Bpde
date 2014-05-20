#ifndef BSOLVER_OMP_H
#define BSOLVER_OMP_H

#include "BArea.h"
#include "BSolver.h"
#include "omp.h"

namespace Bpde {

namespace __bpde_omp {
    static double *H, *Ha;
    static double *x, *y;
    static double *V;
    static double *mu;
    static double *dx_d, *dx_l, *dx_u, *dy_d, *dy_l, *dy_u;
    static double *loc_c, *loc_d, *b;
    static int I, J, T;
    static double dt;
    static int n;
    static double t;
}

class BSolverOmp : public BSolver
{
public:
    BSolverOmp(const BArea& area, int threadsNum);
    virtual ~BSolverOmp();

    double* solve();

    virtual double exec_time();
    int it_num();
private:
    inline void prepareIteration();
    int iterations;
    double time;
    int threadsNum;

    BArea area;
};

}

#endif
