#ifndef BSOLVER_OCL_H
#define BSOLVER_OCL_H

#include "BSolver.h"

#include "CL/cl.hpp"
#include "omp.h"

#include <iostream>
#include "math.h"

namespace Bpde {

class BSolverOcl : public BSolver
{
public:
    BSolverOcl(const BArea& area, std::vector<cl::Device>& devices);
    virtual ~BSolverOcl();

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
    inline cl_int initOpenCL();
    inline cl_int setArgsToExplicitDerivativeKernel();
    inline cl_int setArgsToHydraulicConductivityKernel();

    inline void getMu(double* mu, double H)
    {
        *mu = (H >= area.zc)?area.mu1:area.mu2;
    }

    BArea area;

    double *H, *Ha;
    double *x, *y;
    double *V;
    double *mu;
    double *dx_d, *dx_l, *dx_u, *dy_d, *dy_l, *dy_u;
    double *loc_c, *loc_d, *b, *tmp_v;
    int I, J, T;
    double dt;
    int n;
    double t;

    std::vector<cl::Device> devices;
    cl::Context context;
    cl::CommandQueue commandQueue;
    cl::Program program;
    cl::Kernel explicitDerivativeKernel;
    cl::Kernel hydraulicConductivityKernel;

    cl::Buffer HBuff, HaBuff;
    cl::Buffer xBuff, yBuff;
    cl::Buffer VBuff, muBuff;
    cl::Buffer dx_lBuff, dx_dBuff, dx_uBuff;
    cl::Buffer dy_lBuff, dy_dBuff, dy_uBuff;
    cl::Buffer IBuff, JBuff, TBuff;
    cl::Buffer dtBuff, nBuff, tBuff;
    cl::Buffer zcBuff, zfBuff, kxBuff, kyBuff;
};

} // namespace Bpde

#endif
