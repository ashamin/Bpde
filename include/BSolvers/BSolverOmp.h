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
    static int I, J;
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

    virtual void addExtraIterations(int its);
    virtual void setTimeStep(double dt);

private:    
    inline void prepareIteration();
    inline double Tx(double H)
    {
        return (H >= area.zc)?area.kx*(area.zc - area.zf):((H < area.zf)?0:area.kx*(H - area.zf));
    }

    inline double Ty(double H)
    {
        return (H >= area.zc)?area.ky*(area.zc - area.zf):((H < area.zf)?0:area.ky*(H - area.zf));
    }

    inline void getMu(double* mu, double H)
    {
        *mu = (H >= area.zc)?area.mu1:area.mu2;
    }

    int iterations;
    double time;
    int threadsNum;

    BArea area;


};

}

#endif
