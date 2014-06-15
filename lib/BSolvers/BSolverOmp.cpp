#include "BSolverOmp.h"

#include "stdlib.h"
#include <cstring>

#include <iostream>
#include <fstream>

using namespace std;

namespace Bpde
{

BSolverOmp::BSolverOmp(const BArea& area, int threadsNum)
    :area(area),
     iterations(0),
     time(0),
     threadsNum(threadsNum)
{
    dt = area.dt;

    I = area.I;
    J = area.J;
    n = I*J;

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

    omp_set_dynamic(0);
    omp_set_num_threads(threadsNum);

    t = 0;
    for (int i =0; i<n; i++)
        Ha[i] = b[i] = V[i] = dx_d[i] = dx_l[i] = dx_u[i] = dy_d[i] = dy_l[i] = dy_u[i] =
                mu[i] = loc_c[i] = loc_d[i] = 0;

    tmp_v = new double[n];
}
    
BSolverOmp::~BSolverOmp()
{
    delete[] tmp_v;
}

double* BSolverOmp::solve()
{
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

double BSolverOmp::exec_time()
{
    return time;
}

void BSolverOmp::addExtraIterations(int its)
{
    area.T += its;
}

void BSolverOmp::setTimeStep(double dt)
{
    this->dt = dt;
}

void BSolverOmp::cHydraulicConductivity()
{
    int i, j, k, kT, kH;
    double *dx_l = this->dx_l, *dx_d = this->dx_d, *dx_u = this->dx_u;
    double *dy_l = this->dy_l, *dy_d = this->dy_d, *dy_u = this->dy_u;
    double *H = this->H, *x = this->x, *y = this->y;
    int I = this->I, J = this->J;

#pragma omp parallel for shared(dx_l, dx_d, dx_u) \
firstprivate(I, J) private(j, k) \
schedule(static)
    for (j = 2; j<J-2; j++) {
        k = j*I + 1;
        dx_l[k] = 0;
        dx_d[k] = 1;
        dx_u[k] = -1;

        k = j * I + I - 2;
        dx_l[k] = -1;
        dx_d[k] = 1;
        dx_u[k] = 0;
    }

    // формируем Tx на каждом шаге
#pragma omp parallel for shared(H, dx_l, dx_d, dx_u, x) \
firstprivate(I, J) private(i, j, k) \
schedule(static) collapse(2)
    for (j = 2; j<J-2; j++) {
        for (i = 2; i<I-2; i++){
            k = j*I+i;
            dx_l[k] = Tx((H[k-1] + H[k])/2) /
                    ((x[i] - x[i-1])
                    * ((x[i] + x[i+1])/2 - (x[i] + x[i-1])/2)
                    );
            dx_u[k] = Tx((H[k+1] + H[k])/2) /
                    ((x[i+1] - x[i])
                    * ((x[i] + x[i+1])/2 - (x[i] + x[i-1])/2)
                    );
            dx_d[k] = (-Tx((H[k+1] + H[k])/2) / (x[i+1] - x[i]) -
                        Tx((H[k-1] + H[k])/2) / (x[i] - x[i-1]))
                    / ((x[i] + x[i+1])/2 - (x[i] + x[i-1])/2);
        }
    }

#pragma omp parallel for shared(dy_l, dy_d, dy_u) \
firstprivate(I, J) private(i, kT) \
schedule(static)
    for (i = 2; i < I - 2; i++) {
        int kT = i*J + 1;
        dy_l[kT] = 0;
        dy_d[kT] = 1;
        dy_u[kT] = -1;

        kT = i*J + (J-2);
        dy_l[kT] = -1;
        dy_d[kT] = 1;
        dy_u[kT] = 0;
    }


    // формируем Ty на каждом шаге
#pragma omp parallel for shared(H, dy_l, dy_d, dy_u, y) \
firstprivate(I, J) private(i, j, kT, kH) \
schedule(static) collapse(2)
    for (i = 2; i < I - 2; i++) {
        for (j = 2; j < J - 2; j++) {
            kH = j*I+i;
            kT = i*J + j;
            dy_l[kT] = Ty((H[kH-I] + H[kH])/2) /
                    ((y[j] - y[j-1])
                    * ((y[j] + y[j+1])/2 - (y[j] + y[j-1])/2)
                    );
            dy_u[kT] = Ty((H[kH+I] + H[kH])/2) /
                    ((y[j+1] - y[j])
                    * ((y[j] + y[j+1])/2 - (y[j] + y[j-1])/2)
                    );
            dy_d[kT] =(-Ty((H[kH+I] + H[kH])/2) / (y[j+1] - y[j]) -
                    Ty((H[kH-I] + H[kH])/2) / (y[j] - y[j-1]))
                    / ((y[j] + y[j+1])/2 - (y[j] + y[j-1])/2);
        }
    }
}

void BSolverOmp::cSupporting()
{
    int i, j;
    double *V = this->V, *x = this->x, *y = this->y;
    double *mu = this->mu, *H = this->H;
    int I = this->I, J = this->J, n = this->n;
    double t = this->t;
    // пересчитываем функцию V а каждом шаге
#pragma omp parallel for shared(V) \
firstprivate(I, J, t) private(i, j) \
schedule(static) collapse(2)
    for (j = 1; j<J-1; j++)
        for (i = 1; i<I-1; i++)
            V[i + I*j] = area.V(x[i], y[j], t);

    // вычисляем зачения mu на всей области (с буфером)
#pragma omp parallel for shared(mu) \
firstprivate(I, n) private(i) \
schedule(static)
    for (i = I; i<n-I; i++)
        getMu(&mu[i], H[i]);
}

void BSolverOmp::cExplicitDerivative()
{

    int i, j, kT, kH;
    double *H = this->H, *Ha = this->Ha, *V = this->V, *mu = this->mu;
    double *dx_l = this->dx_l, *dx_d = this->dx_d, *dx_u = this->dx_u;
    double *dy_l = this->dy_l, *dy_d = this->dy_d, *dy_u = this->dy_u;
    int I = this->I, J = this->J;

#pragma omp parallel for shared(H, dx_l, dx_d, dx_u, dy_l, dy_d, dy_u, V, mu) \
firstprivate(I, J) private(i, j, kT, kH) \
schedule(static) collapse(2)
    for (i = 1; i < I - 1; i++)
        for (j = 1; j < J - 1; j++) {
            kT = i*J + j;
            kH = j*I + i;
            Ha[kH] =(
                     dx_l[kH]*H[kH-1] + dx_d[kH]*H[kH] + dx_u[kH]*H[kH+1] +
                     dy_l[kT]*H[kH-I] + dy_d[kT]*H[kH] + dy_u[kT]*H[kH+I] +
                     V[kH]
                    )
                    / mu[kH];
        }
}

void BSolverOmp::cImplicitTDMAs()
{
    int i, j, k, kT, kH;
    double tmp;
    double *H = this->H, *Ha = this->Ha, *mu = this->mu, *b = this->b;
    double *dx_l = this->dx_l, *dx_d = this->dx_d, *dx_u = this->dx_u;
    double *dy_l = this->dy_l, *dy_d = this->dy_d, *dy_u = this->dy_u;
    double *loc_c = this->loc_c, *loc_d = this->loc_d;
    int I = this->I, J = this->J;
    double dt = this->dt;
    double *tmp_v = this->tmp_v;

    // неявная прогонка по X
#pragma omp parallel for shared(Ha, b, mu, dx_d) \
firstprivate(I, J, dt) private(i, j, k, tmp) \
schedule(static) collapse(2)
    for (i = 2; i < I - 2; i++)
        for (j = 2; j < J - 2; j++) {
            k = j*I + i;
            tmp = mu[k] / dt;
            dx_d[k] -= tmp;
            b[k] = - Ha[k] * tmp;
        }

#pragma omp parallel for shared(dx_l, dx_d, dx_u, tmp_v, b, loc_c, loc_d) \
firstprivate(I, J) private(j, k) \
schedule(static)
    for (j = 2; j<J-2; j++){
        k = j*I+1;
        TDMA(&dx_l[k], &dx_d[k], &dx_u[k], &tmp_v[k], &b[k], I-2, 1, &loc_c[k], &loc_d[k]);
    }

    // неявная прогонка по Y
#pragma omp parallel for shared(tmp_v, mu, dy_d) \
firstprivate(I, J, dt) private(i, j, kT, kH, tmp) \
schedule(static) collapse(2)
    for (i = 2; i < I - 2; i++)
        for (j = 2; j < J - 2; j++) {
            kT = i*J + j;
            kH = j*I + i;
            tmp = mu[kH] / dt;
            dy_d[kT] -= tmp;
            tmp_v[kH] = -tmp_v[kH] * tmp;
    }

#pragma omp parallel for shared(tmp_v) \
firstprivate(I, J) private(i, kH) \
schedule(static)
    for (i = 2; i< I - 2; i++) {
        int kH = I + i;
        tmp_v[kH] = 0;
        kH = (J-2)*I + i;
        tmp_v[kH] =  0;
    }

#pragma omp parallel for shared(dy_l, dy_d, dy_u, Ha, tmp_v, loc_c, loc_d) \
firstprivate(I, J) private(i, kH, kT) \
schedule(static)
    for (i = 2; i < I - 2; i++) {
        kH = I + i;
        kT = i*J+1;
        TDMA_t(&dy_l[kT], &dy_d[kT], &dy_u[kT], &Ha[kH], &tmp_v[kH], J-2, I, &loc_c[kT], &loc_d[kT]);
    }
}

void BSolverOmp::cNextLayer()
{
    int i, j;
    double *H = this->H, *Ha = this->Ha;
    int I = this->I, J = this->J;
    double dt = this->dt;
#pragma omp parallel for shared(H, Ha) \
firstprivate(I, J) private(i, j) \
schedule(static)
    for (i = 2; i < I - 2; i++)
        for (j = 2; j < J - 2; j++)
        H[j*I + i] = H[j*I + i] + dt*Ha[j*I + i];
}

} // namespace Bpde
