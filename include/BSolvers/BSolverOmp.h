#ifndef BSOLVER_OMP_H
#define BSOLVER_OMP_H

#include "omp.h"

#include "BSolver.h"
#include "BLogger/BLogger.h"

namespace Bpde {

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

protected:
    inline virtual void cHydraulicConductivity();
    inline virtual void cSupporting();
    inline virtual void cExplicitDerivative();
    inline virtual void cImplicitTDMAs();
    inline virtual void cNextLayer();

private:    
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

    double *tmp_v;

    BArea area;
};

}

#endif
