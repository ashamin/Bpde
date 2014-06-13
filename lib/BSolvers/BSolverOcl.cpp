#include "BSolverOcl.h"

#include "stdlib.h"
#include <cstring>

#include <iostream>
#include <fstream>

#define PTR_FLAG CL_MEM_USE_HOST_PTR

using namespace std;

namespace Bpde
{

BSolverOcl::BSolverOcl(const BArea& area, std::vector<cl::Device>& devices)
    :area(area),
     devices(devices), context(), commandQueue(), program(),
     explicitDerivativeKernel(), hydraulicConductivityKernel(),
     t(0), dt(area.dt), I(area.I), J(area.J), n(area.I * area.J),
     H(NULL), Ha(NULL), b(NULL), V(NULL), mu(NULL), loc_c(NULL), loc_d(NULL),
     dx_d(NULL), dx_l(NULL), dx_u(NULL), dy_d(NULL), dy_l(NULL), dy_u(NULL)
{
    H     = area.H;
    x     = area.x;
    y     = area.y;

    Ha    = new double[n];
    b     = new double[n];
    V     = new double[n];
    dx_d  = new double[n];
    dx_l  = new double[n];
    dx_u  = new double[n];
    dy_d  = new double[n];
    dy_l  = new double[n];
    dy_u  = new double[n];
    mu    = new double[n];

    loc_c = new double[n];
    loc_d = new double[n];

    for (int i =0; i<n; i++)
        Ha[i] = b[i] = V[i] = dx_d[i] = dx_l[i] = dx_u[i] = dy_d[i] = dy_l[i] = dy_u[i] =
                mu[i] = loc_c[i] = loc_d[i] = 0;
    tmp_v = new double[n];
}

BSolverOcl::~BSolverOcl()
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

cl_int BSolverOcl::initOpenCL()
{
    cl_int ret = 0, err = 0;
    context = cl::Context(devices);
    try {
//        std::ifstream file("kernels.cl");
        std::ifstream file("/home/ashamin/diploma/src/Bpde/lib/BSolvers/kernels.cl");
        std::string programString(std::istreambuf_iterator<char>(file),
            (std::istreambuf_iterator<char>()));
        cl::Program::Sources source(1, std::make_pair(programString.c_str(),
            programString.length()+1));

        program = cl::Program(context, source, &ret);
        err |= ret;
        err |= program.build(devices);

        commandQueue = cl::CommandQueue(context, 0, &ret);
        err |= ret;

        explicitDerivativeKernel = cl::Kernel(program, "explicitDerivative", &ret);
        err |= ret;
        hydraulicConductivityKernel = cl::Kernel(program, "hydraulicConductivityKernel", &ret);
        err |= ret;
    } catch(...)
    {
        std::cout << "exception were threwed" << std::endl;
        return -1;
    }

    return err;
}

double* BSolverOcl::solve()
{
    if (initOpenCL() < 0)
    {
        std::cout << "OpenCL does not initialized" << std::endl;
        exit(-1);
    }

    cl_int err;

    double* tmp_v = new double[n];

    time = omp_get_wtime();

    while (iterations < area.T){
        iterate();
        t += dt;
        iterations++;
        if (iterations % 10000 == 0)
            std::cout << iterations << std::endl;
    }

    time = omp_get_wtime() - time;

    return H;
}

double BSolverOcl::exec_time()
{
    return time;
}

cl_int BSolverOcl::setArgsToExplicitDerivativeKernel()
{
    cl_int err = 0, ret = 0;
    HBuff = cl::Buffer(context,
            CL_MEM_READ_ONLY | PTR_FLAG,
            sizeof(double)*n, H, &ret);
    err |= ret;
    VBuff = cl::Buffer(context,
            CL_MEM_READ_ONLY | PTR_FLAG,
            sizeof(double)*n, V, &ret);
    err |= ret;
    muBuff = cl::Buffer(context,
            CL_MEM_READ_ONLY | PTR_FLAG,
            sizeof(double)*n, mu, &ret);
    err |= ret;
    dx_lBuff = cl::Buffer(context,
            CL_MEM_READ_ONLY | PTR_FLAG,
            sizeof(double)*n, dx_l, &ret);
    err |= ret;
    dx_dBuff = cl::Buffer(context,
            CL_MEM_READ_ONLY | PTR_FLAG,
            sizeof(double)*n, dx_d, &ret);
    err |= ret;
    dx_uBuff = cl::Buffer(context,
            CL_MEM_READ_ONLY | PTR_FLAG,
            sizeof(double)*n, dx_u, &ret);
    err |= ret;
    dy_lBuff = cl::Buffer(context,
            CL_MEM_READ_ONLY | PTR_FLAG,
            sizeof(double)*n, dy_l, &ret);
    err |= ret;
    dy_dBuff = cl::Buffer(context,
            CL_MEM_READ_ONLY | PTR_FLAG,
            sizeof(double)*n, dy_d, &ret);
    err |= ret;
    dy_uBuff = cl::Buffer(context,
            CL_MEM_READ_ONLY | PTR_FLAG,
            sizeof(double)*n, dy_u, &ret);
    err |= ret;
    IBuff = cl::Buffer(context,
            CL_MEM_READ_ONLY | PTR_FLAG,
            sizeof(int), &I, &ret);
    err |= ret;
    JBuff = cl::Buffer(context,
            CL_MEM_READ_ONLY | PTR_FLAG,
            sizeof(int), &J, &ret);
    err |= ret;
    HaBuff = cl::Buffer(context,
            CL_MEM_WRITE_ONLY | PTR_FLAG,
            sizeof(double)*n, Ha, &ret);
    err |= ret;

    err |= explicitDerivativeKernel.setArg(0, HBuff);
    err |= explicitDerivativeKernel.setArg(1, VBuff);
    err |= explicitDerivativeKernel.setArg(2, muBuff);
    err |= explicitDerivativeKernel.setArg(3, dx_lBuff);
    err |= explicitDerivativeKernel.setArg(4, dx_dBuff);
    err |= explicitDerivativeKernel.setArg(5, dx_uBuff);
    err |= explicitDerivativeKernel.setArg(6, dy_lBuff);
    err |= explicitDerivativeKernel.setArg(7, dy_dBuff);
    err |= explicitDerivativeKernel.setArg(8, dy_uBuff);
    err |= explicitDerivativeKernel.setArg(9, IBuff);
    err |= explicitDerivativeKernel.setArg(10, JBuff);
    err |= explicitDerivativeKernel.setArg(11, HaBuff);
}

cl_int BSolverOcl::setArgsToHydraulicConductivityKernel()
{
    cl_int err = 0, ret = 0;
    HBuff = cl::Buffer(context,
            CL_MEM_READ_ONLY | PTR_FLAG,
            sizeof(double)*n, H, &ret);
    err |= ret;
    xBuff = cl::Buffer(context,
            CL_MEM_READ_ONLY | PTR_FLAG,
            sizeof(double)*I, x, &ret);
    err |= ret;
    yBuff = cl::Buffer(context,
            CL_MEM_READ_ONLY | PTR_FLAG,
            sizeof(double)*J, y, &ret);
    err |= ret;
    IBuff = cl::Buffer(context,
            CL_MEM_READ_ONLY | PTR_FLAG,
            sizeof(int), &I, &ret);
    err |= ret;
    JBuff = cl::Buffer(context,
            CL_MEM_READ_ONLY | PTR_FLAG,
            sizeof(int), &J, &ret);
    err |= ret;
    zcBuff = cl::Buffer(context,
            CL_MEM_READ_ONLY | PTR_FLAG,
            sizeof(double), &area.zc, &ret);
    err |= ret;
    zfBuff = cl::Buffer(context,
            CL_MEM_READ_ONLY | PTR_FLAG,
            sizeof(double), &area.zf, &ret);
    err |= ret;
    kxBuff = cl::Buffer(context,
            CL_MEM_READ_ONLY | PTR_FLAG,
            sizeof(double), &area.kx, &ret);
    err |= ret;
    kyBuff = cl::Buffer(context,
            CL_MEM_READ_ONLY | PTR_FLAG,
            sizeof(double), &area.ky, &ret);
    err |= ret;
    dx_lBuff = cl::Buffer(context,
            CL_MEM_WRITE_ONLY | PTR_FLAG,
            sizeof(double)*n, dx_l, &ret);
    err |= ret;
    dx_dBuff = cl::Buffer(context,
            CL_MEM_WRITE_ONLY | PTR_FLAG,
            sizeof(double)*n, dx_d, &ret);
    err |= ret;
    dx_uBuff = cl::Buffer(context,
            CL_MEM_WRITE_ONLY | PTR_FLAG,
            sizeof(double)*n, dx_u, &ret);
    err |= ret;
    dy_lBuff = cl::Buffer(context,
            CL_MEM_WRITE_ONLY | PTR_FLAG,
            sizeof(double)*n, dy_l, &ret);
    err |= ret;
    dy_dBuff = cl::Buffer(context,
            CL_MEM_WRITE_ONLY | PTR_FLAG,
            sizeof(double)*n, dy_d, &ret);
    err |= ret;
    dy_uBuff = cl::Buffer(context,
            CL_MEM_WRITE_ONLY | PTR_FLAG,
            sizeof(double)*n, dy_u, &ret);
    err |= ret;

    err |= hydraulicConductivityKernel.setArg(0, HBuff);
    err |= hydraulicConductivityKernel.setArg(1, xBuff);
    err |= hydraulicConductivityKernel.setArg(2, yBuff);
    err |= hydraulicConductivityKernel.setArg(3, IBuff);
    err |= hydraulicConductivityKernel.setArg(4, JBuff);
    err |= hydraulicConductivityKernel.setArg(5, zcBuff);
    err |= hydraulicConductivityKernel.setArg(6, zfBuff);
    err |= hydraulicConductivityKernel.setArg(7, kxBuff);
    err |= hydraulicConductivityKernel.setArg(8, kyBuff);
    err |= hydraulicConductivityKernel.setArg(9, dx_lBuff);
    err |= hydraulicConductivityKernel.setArg(10, dx_dBuff);
    err |= hydraulicConductivityKernel.setArg(11, dx_uBuff);
    err |= hydraulicConductivityKernel.setArg(12, dy_lBuff);
    err |= hydraulicConductivityKernel.setArg(13, dy_dBuff);
    err |= hydraulicConductivityKernel.setArg(14, dy_uBuff);
}

void BSolverOcl::addExtraIterations(int its)
{
    area.T += its;
}

void BSolverOcl::setTimeStep(double dt)
{
    this->dt = dt;
}

void BSolverOcl::cHydraulicConductivity()
{
    cl_int err = 0;
    int k = 0;
    // формируем Tx и Ty на каждом шаге
    setArgsToHydraulicConductivityKernel();

    err = commandQueue.enqueueNDRangeKernel(hydraulicConductivityKernel,
                cl::NDRange(2, 2), cl::NDRange(J-4, I-4));

    err = commandQueue.enqueueReadBuffer(dx_lBuff, CL_TRUE, 0,
                sizeof(double)*n, static_cast<void*>(dx_l));
    err = commandQueue.enqueueReadBuffer(dx_dBuff, CL_TRUE, 0,
                sizeof(double)*n, static_cast<void*>(dx_d));
    err = commandQueue.enqueueReadBuffer(dx_uBuff, CL_TRUE, 0,
                sizeof(double)*n, static_cast<void*>(dx_u));

    err = commandQueue.enqueueReadBuffer(dy_lBuff, CL_TRUE, 0,
                sizeof(double)*n, static_cast<void*>(dy_l));
    err = commandQueue.enqueueReadBuffer(dy_dBuff, CL_TRUE, 0,
                sizeof(double)*n, static_cast<void*>(dy_d));
    err = commandQueue.enqueueReadBuffer(dy_uBuff, CL_TRUE, 0,
                sizeof(double)*n, static_cast<void*>(dy_u));

    for (int j = 2; j<J-2; j++) {
        k = j*I + 1;
        dx_l[k] = 0;
        dx_d[k] = 1;
        dx_u[k] = -1;

        k = j * I + I - 2;
        dx_l[k] = -1;
        dx_d[k] = 1;
        dx_u[k] = 0;
    }

    for (int i = 2; i < I - 2; i++) {
        int kT = i*J + 1;
        dy_l[kT] = 0;
        dy_d[kT] = 1;
        dy_u[kT] = -1;

        kT = i*J + (J-2);
        dy_l[kT] = -1;
        dy_d[kT] = 1;
        dy_u[kT] = 0;
    }
}

void BSolverOcl::cSupporting()
{
    // пересчитываем функцию V а каждом шаге
    for (int j = 1; j<J-1; j++)
        for (int i = 1; i<I-1; i++)
            V[i + I*j] = area.V(x[i], y[j], t);

    // вычисляем зачения mu на всей области (с буфером)
    for (int i = I; i<n-I; i++)
        getMu(&mu[i], H[i]);
}

void BSolverOcl::cExplicitDerivative()
{
    cl_int err = 0;
    setArgsToExplicitDerivativeKernel();

    err = commandQueue.enqueueNDRangeKernel(explicitDerivativeKernel,
                cl::NDRange(1, 1), cl::NDRange(I-2, J-2));

    err = commandQueue.enqueueReadBuffer(HaBuff, CL_TRUE, 0,
                                         sizeof(double)*n, static_cast<void*>(Ha));
}

void BSolverOcl::cImplicitTDMAs()
{
    // неявная прогонка по X
    for (int i = 2; i < I - 2; i++)
        for (int j = 2; j < J - 2; j++) {
            int k = j*I + i;
            double tmp = mu[k] / dt;
            dx_d[k] -= tmp;
            b[k] = - Ha[k] * tmp;
        }

    for (int j = 2; j<J-2; j++){
        int k = j*I+1;
        TDMA(&dx_l[k], &dx_d[k], &dx_u[k], &tmp_v[k], &b[k], I-2, 1, &loc_c[k], &loc_d[k]);
    }

    // неявная прогонка по Y
    for (int i = 2; i < I - 2; i++)
        for (int j = 2; j < J - 2; j++) {
            int kT = i*J + j;
            int kH = j*I + i;
            double tmp = mu[kH] / dt;
            dy_d[kT] -= tmp;
            tmp_v[kH] = -tmp_v[kH] * tmp;
    }

    for (int i = 2; i< I - 2; i++) {
        int kH = I + i;
        tmp_v[kH] = 0;
        kH = (J-2)*I + i;
        tmp_v[kH] =  0;
    }

    for (int i = 2; i < I - 2; i++) {
        int kH = I + i;
        int kT = i*J+1;
        TDMA_t(&dy_l[kT], &dy_d[kT], &dy_u[kT], &Ha[kH], &tmp_v[kH], J-2, I, &loc_c[kT], &loc_d[kT]);
    }
}

void BSolverOcl::cNextLayer()
{
    for (int i = 2; i < I - 2; i++)
        for (int j = 2; j < J - 2; j++)
            H[j*I + i] = H[j*I + i] + dt*Ha[j*I + i];
}


} // namespace Bpde
