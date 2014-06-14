#ifndef BSOLVER_H
#define BSOLVER_H

#include <iostream>
#include <cmath>

#include "BArea.h"

namespace Bpde {

namespace ParallelizationMethod
{
enum ParallelizationMethod
{
    NONE = 0,
    OPENMP = 1,
    OPENCL = 2
};
} // namespace ParallelizationMethod

/**
 *
 */
class BSolver
{
public:
    BSolver()
        :iterations(0), threadsNum(0), time(0), t(0), n(0), dt(0), I(0), J(0),
          H(NULL), Ha(NULL), x(NULL), y(NULL), V(NULL), mu(NULL), b(NULL),
          dx_d(NULL), dx_l(NULL), dx_u(NULL), dy_d(NULL), dy_l(NULL), dy_u(NULL),
          loc_c(NULL), loc_d(NULL)
    {}
    virtual ~BSolver()
    {
        delete[] Ha;
        delete[] b;
        delete[] V;
        delete[] dx_d;
        delete[] dx_l;
        delete[] dx_u;
        delete[] dy_d;
        delete[] dy_l;
        delete[] dy_u;
        delete[] mu;
        delete[] loc_c;
        delete[] loc_d;
    }
    virtual double* solve(){}
    inline virtual void iterate() {
        cHydraulicConductivity();
        cSupporting();
        cExplicitDerivative();
        cImplicitTDMAs();
        cNextLayer();
    }
    virtual double exec_time(){}
    virtual void addExtraIterations(int its){}
    virtual void setTimeStep(double dt){}
protected:
    inline virtual void cHydraulicConductivity() {}
    inline virtual void cSupporting() {}
    inline virtual void cExplicitDerivative() {}
    inline virtual void cImplicitTDMAs() {}
    inline virtual void cNextLayer() {}

    double *H, *Ha;
    double *x, *y;
    double *V;
    double *mu;
    double *dx_d, *dx_l, *dx_u, *dy_d, *dy_l, *dy_u;
    double *loc_c, *loc_d, *b;
    int I, J;
    double dt;
    int n;
    double t;

    int iterations;
    double time;
    int threadsNum;
};

void TDMA(const double* a, const double* b, const double* c,
        double* x, const double* d, int n, int step, double *loc_c, double *loc_d);

/**
 * @brief TDMA_t
 * @param a
 * @param b
 * @param c
 * @param x
 * @param d
 * @param n
 * @param step
 * @param loc_c
 * @param loc_d
 *
 * Трехдиагональная матрица ориентирована на транспонированную матрицу в отличие от TDMA,
 *  где трехдиагональная матрица ориентирована на исходную матрицу.
 */
void TDMA_t(const double* a, const double* b, const double* c,
        double* x, const double* d, int n, int step, double *loc_c, double *loc_d);

} // namespace Bpde

#endif
