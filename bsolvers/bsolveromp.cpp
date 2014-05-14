#include "bsolveromp.h"

#include "stdlib.h"
#include <cstring>

#include <iostream>
#include <fstream>

using namespace std;

namespace Bpde
{

BSolverOmp::BSolverOmp(BArea* area)
{
    using namespace __bpde_omp;

    this->area = area;
    t = 0;
    dt = area->dt;

    I = area->I;
    J = area->J;
    n = I*J;

    H     = area->H;
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

    iterations = 0;
}
    
BSolverOmp::~BSolverOmp()
{
    using namespace __bpde_omp;
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

void BSolverOmp::prepareIteration()
{
    using namespace __bpde_omp;
    // пересчитываем функцию V а каждом шаге
    for (int j = 1; j<J-1; j++)
        for (int i = 1; i<I-1; i++)
            V[i + I*j] = area->V(area->x[i], area->y[j], t);

    // формируем Tx на каждом шаге
    for (int j = 2; j<J-2; j++) {
        int k = j*I + 1;
        dx_l[k] = 0;
        dx_d[k] = 1;
        dx_u[k] = -1;

        for (k = j*I+2; k<j*I+I-2; k++) {
            dx_l[k] = Tx((H[k-1] + H[k])/2) / (area->hx * area->hx);
            dx_u[k] = Tx((H[k+1] + H[k])/2) / (area->hx * area->hx);
            dx_d[k] = (-Tx((H[k+1] + H[k])/2) /
                    area->hx - Tx((H[k-1] + H[k])/2) / area->hx) / area->hx;
        }

        k = j * I + I - 2;
        dx_l[k] = -1;
        dx_d[k] = 1;
        dx_u[k] = 0;
    }

    // формируем Ty на каждом шаге
    for (int i = 2; i < I - 2; i++) {
        int kT = i*J + 1;
        dy_l[kT] = 0;
        dy_d[kT] = 1;
        dy_u[kT] = -1;

        for (int j = 2; j < J - 2; j++) {
            int kH = j*I+i;
            kT = i*J + j;
            dy_l[kT] = Ty((H[kH-I] + H[kH])/2) / (area->hy * area->hy);
            dy_u[kT] = Ty((H[kH+I] + H[kH])/2) / (area->hy * area->hy);
            dy_d[kT] =(-Ty((H[kH+I] + H[kH])/2) /
                    area->hy - Ty((H[kH-I] + H[kH])/2) / area->hy) / area->hy;
        }

        kT = i*J + (J-2);
        dy_l[kT] = -1;
        dy_d[kT] = 1;
        dy_u[kT] = 0;
    }
    

    // вычисляем зачения mu на всей области (с буфером)
    for (int i = I; i<n-I; i++)
        getMu(&mu[i], H[i]);
}

double* BSolverOmp::solve()
{
    using namespace __bpde_omp;

    double* tmp_v = new double[n];

    while (t<(area->dt*(area->T-1))){

        prepareIteration();

        int k = 0;
        for (int i = 1; i < I - 1; i++)
            for (int j = 1; j < J - 1; j++) {
                int kT = i*J + j;
                int kH = j*I + i;
                Ha[kH] =(
                         dx_l[kH]*H[kH-1] + dx_d[kH]*H[kH] + dx_u[kH]*H[kH+1] +
                         dy_l[kT]*H[kH-I] + dy_d[kT]*H[kH] + dy_u[kT]*H[kH+I] +
                         V[kH]
                        )
                        / mu[kH];
         }

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

        for (int i = 2; i< I -2; i++) {
            int kH = I + i;
            tmp_v[kH] = 0;
            kH = (J-2)*I + i;
            tmp_v[kH] =  0;
        }

        // исправить это!
        for (int i = 2; i < I - 2; i++) {
            int kH = I + i;
            int kT = i*J+1;
            TDMA_t(&dy_l[kT], &dy_d[kT], &dy_u[kT], &Ha[kH], &tmp_v[kH], J-2, I, &loc_c[kT], &loc_d[kT]);
        }

        for (int i = 2; i < I - 2; i++)
            for (int j = 2; j < J - 2; j++)
            H[j*I + i] = H[j*I + i] + dt*Ha[j*I + i];


        t += dt;
        iterations++;
        if (iterations % 10000 == 0)
            std::cout << iterations << std::endl;
    }

    log_matrix("H", Ha, I, J);
    std::cout << std::endl << std::endl;
    log_matrix("H", H, I, J);

    return H;
}

double BSolverOmp::exec_time()
{
    return 0.0;
}

} // namespace Bpde
