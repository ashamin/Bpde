#include "bsolveromp.h"

#include "stdlib.h"
#include <cstring>

using namespace std;

namespace Boussinesq
{

BSolverOmp::BSolverOmp(BArea* area, const double epsilon, const int maxit)
{


    this->area = area;
    this->epsilon = epsilon;
    this->maxit = maxit;
    t = 0;
    dt = area->dt;

    I = area->I + 2;
    J = area->J + 2;
    n = I*J;

    x = new double[area->I];
    y = new double[area->J];

    x[0] = y[0] = 0;
    for (int i = 1; i<area->I; i++)
        x[i] = x[i-1] + area->hx;
    for (int j = 1; j<area->J; j++)
        y[j] = y[j-1] + area->hy;


    H   = new double[n];
    Ha  = new double[n];

//    memset(H, 0, n*sizeof(double));
//    memset(Ha, 0, n*sizeof(double));

//    for (int j = 1; j<J-1; j++)
//        for (int i = 1; i<I-1; i++)
//            H[i + I*j] = area->answer(x[i-1], y[j-1], 0);
//    for (int i = 0; i<n; i++)
//        H[i] = area->H[i];

    H = area->H;

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

    memset(b, 0, n*sizeof(double));
    memset(V, 0, n*sizeof(double));
    memset(dx_d, 0, n*sizeof(double));
    memset(dx_l, 0, n*sizeof(double));
    memset(dx_u, 0, n*sizeof(double));
    memset(dy_d, 0, n*sizeof(double));
    memset(dy_l, 0, n*sizeof(double));
    memset(dy_u, 0, n*sizeof(double));
    memset(mu, 0, n*sizeof(double));
    memset(loc_c, 0, n*sizeof(double));
    memset(loc_d, 0, n*sizeof(double));

    iterations = 0;
}
    
BSolverOmp::~BSolverOmp()
{

}

double BSolverOmp::borderValue(double x, double y, double t)
{
    return area->answer(x, y, t);
}

void BSolverOmp::recomputeBorderValues()
{

}

void BSolverOmp::prepareIteration()
{
    // пересчитываем функцию V а каждом шаге
    for (int j = 1; j<J-1; j++)
        for (int i = 1; i<I-1; i++)
            V[i + I*j] = area->V(x[i-1], y[j-1], t);

    // формируем Tx на каждом шаге
    for (int j = 2; j<J-2; j++) {
        int k = j*I + 1;
        dx_l[k] = 0;
        dx_d[k] = 1;
        dx_u[k] = -1;

        for (k = j*I+2; k<j*I+I-2; k++) {
//            if (k == 126)
//                std::cout << "olol" << std::endl;
//            dx_l[k]  = Tx(H[k]) / (area->hx * area->hx) -
//                    (Tx((H[k+1] + H[k])/2) - Tx((H[k-1] + H[k])/2)) / (2 * area->hx);
//            dx_u[k] = Tx(H[k]) / (area->hx * area->hx) +
//                    (Tx((H[k+1] + H[k])/2) - Tx((H[k-1] + H[k])/2)) / (2 * area->hx);
//            dx_d[k] = -2 * Tx(H[k]) / (area->hx * area->hx);
            dx_l[k] = Tx((H[k-1] + H[k])/2) / (area->hx * area->hx);
            dx_u[k] = Tx((H[k+1] + H[k])/2) / (area->hx * area->hx);
            dx_d[k] = (-Tx((H[k+1] + H[k])/2) / area->hx - Tx((H[k-1] + H[k])/2) / area->hx) / area->hx;
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
//            if (kT == 98)
//                std::cout << "olol" << std::endl;
//            dy_l[kT] = Ty(H[kH]) / (area->hy * area->hy) -
//                    (Ty((H[kH+I] + H[kH])/2) - Ty((H[kH-I] + H[kH])/2)) / (2 * area->hy);
//            dy_u[kT] = Ty(H[kH]) / (area->hy * area->hy) +
//                    (Ty((H[kH+I] + H[kH])/2) - Ty((H[kH-I] + H[kH])/2)) / (2 * area->hy);
//            dy_d[kT] = -2 * Tx(H[kH]) / (area->hy * area->hy);
            dy_l[kT] = Ty((H[kH-I] + H[kH])/2) / (area->hy * area->hy);
            dy_u[kT] = Ty((H[kH+I] + H[kH])/2) / (area->hy * area->hy);
            dy_d[kT] =(-Ty((H[kH+I] + H[kH])/2) / area->hy - Ty((H[kH-I] + H[kH])/2) / area->hy) / area->hy;
        }

        kT = i*J + (J-2);
        dy_l[kT] = -1;
        dy_d[kT] = 1;
        dy_u[kT] = 0;
    }
    

    // вычисляем зачения mu на всей области (с буфером)
    for (int i = I; i<n-I; i++)
        getMu(&mu[i], H[i]);

    // обнуляем элементы mu, находящиеся в области относящейся к буферу
//    for (int i = I; i<n; i+=I)
//        mu[i] = 0;
//    for (int i = I-1; i<=n; i+=I)
//        mu[i] = 0;

}

double* BSolverOmp::solve()
{
  
    double* tmp_v = new double[n];
    int s         = (int)sqrt(n);

//    int tk = 17*I+324;
    int tk = 7*I + 12;
//    std::cout << std::endl << area->hx << std::endl << std::endl;
//    std::cout << std::endl << H[tk] << "\t";
//    log_matrix("H", H, I, J);


    while (t<(area->dt*(area->T-1))){

        prepareIteration();

//        if (iterations % 100 == 0)
//            std::cout << iterations << "\t" << H[tk] << "\t";

//        log_matrix("H", H, I, J);
//        log_diags_as_3dmatrix("DX", dx_l, dx_d, dx_u, n, I);
//        log_diags_as_3dmatrix("DY", dy_l, dy_d, dy_u, n, J);
//        log_vector("MU", mu, n);
//        log_matrix("V", V, I, J);
//        log_matrix("HA", Ha, I, J);


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

//                if (kH == 126) {
//                    std::cout << "--------------------------------------------" << std::endl;
//                    std::cout << dx_l[kH] << "\t" << dx_d[kH] << "\t" << dx_u[kH] << std::endl;
//                    std::cout << dy_l[kT] << "\t" << dy_d[kT] << "\t" << dy_u[kT] << std::endl;
//                    std::cout << V[kH] << std::endl;
//                    std::cout << mu[kH] << std::endl;
//                    std::cout << "--------------------------------------------" << std::endl;
//                }

//                if (H[kH] > 180){
//                    std::cout << "AAAAAAAAAAAaa ------- " << kH << std::endl;
//                    exit(0);
//                }

                if (kH == 437)
                    std::cout << Ha[kH] << std::endl;
         }

//        log_matrix("HA", Ha, I, J);

        // неявная прогонка по X
//        for (int i = I; i<n-I; i++) {
        for (int i = 2; i < I - 2; i++)
            for (int j = 2; j < J - 2; j++) {
                int k = j*I + i;
                double tmp = mu[k] / dt;
                dx_d[k] -= tmp;
                b[k] = - Ha[k] * tmp;
            }


//        log_diags_as_3dmatrix("DX", dx_l, dx_d, dx_u, n, I);
//        log_vector("B", b, n);

        for (int j = 2; j<J-2; j++){
            int k = j*I+1;
            TDMA(&dx_l[k], &dx_d[k], &dx_u[k], &tmp_v[k], &b[k], I-2, 1, &loc_c[k], &loc_d[k]);
        }

//        log_matrix("TMP", tmp_v, I, J);

        // неявная прогонка по Y
        // возмонжо адресация mu должна быть по kT,
        //  а текущая ошибка по такой адресации из-за того, что mu заполнено не полностью
        // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! -- скорее всего ошибка тут!!!!
        for (int i = 2; i < I - 2; i++)
            for (int j = 2; j < J - 2; j++) {
                int kT = i*J + j;
                int kH = j*I + i;
                double tmp = mu[kH] / dt;
                dy_d[kT] -= tmp;
                tmp_v[kH] = -tmp_v[kH] * tmp;
        }

//        log_diags_as_3dmatrix("DY", dy_l, dy_d, dy_u, n, J);

//        log_vector("TMP_V", tmp_v, n);

//        log_matrix("HA_TDMA", Ha, I, J);

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

//        log_vector("TMP_V", tmp_v, n);
//        log_matrix("HA_TDMA", Ha, I, J);

        // кроме границ
//        for (int i = I; i<n-I; i++)
//            H[i] = H[i]+dt*Ha[i];
        // не учитываем значения на границах. ГРАНИЦЫ КОНСТАНТНЫЕ! ИСПРАВИТЬ!

//        std::cout << iterations << "\t" <<H[126] << "\t" << Ha[126] << std::endl;

        for (int i = 2; i < I - 2; i++)
            for (int j = 2; j < J - 2; j++)
            H[j*I + i] = H[j*I + i] + dt*Ha[j*I + i];


//        if (iterations % 100 == 0)
//            std::cout << H[tk] << std::endl << std::endl;

//        log_matrix("H", H, I, J);

        t += dt;
        iterations++;
        if (iterations % 100 == 0)
//            break;
            std::cout << "IIIITeration --- " << iterations << std::endl;


//        break;

    }


//    std::cout << I << "\t" << J << std::endl;
    log_matrix("Ha", H, I, J);

//    std::cout << iterations << std::endl;

    return H;
}

} // namespace Boussinesq
