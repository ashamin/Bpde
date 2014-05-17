#include "BSolver.h"


namespace Bpde
{
//// к чему это???? wtf
/////!!!!!! Исправить и отладить!!!!!!
/// 2)  loc_d[i] = (d[i] - loc_d[i-step]*a[i]) / tmp;----->>>    loc_d[i] = (d[i*step] - loc_d[i-step]*a[i]) / tmp;
void TDMA(const double *a, const double *b, const double *c,
        double *x, const double *d, int n, int step, double *loc_c, double *loc_d)
{
    double tmp = 0;
    int i = 0;

    loc_c[0] = c[0] / b[0];  // c[0] - всегда 1  (-1)есть
    loc_d[0] = d[0] / b[0]; // d[0] - всегда 0 !!!

    for (i = step; i<n*step; i+=step){
        tmp = b[i] - loc_c[i-step]*a[i];
        loc_c[i] = c[i] / tmp;
        loc_d[i] = (d[i] - loc_d[i-step]*a[i]) / tmp;
    }

    x[n*step-step] = loc_d[n*step-step];
    for (i = (n-2)*step; i>=0; i-=step)
        x[i] = loc_d[i] - loc_c[i]*x[i+step];

    return;
}

void TDMA_t(const double *a, const double *b, const double *c,
        double *x, const double *d, int n, int step, double *loc_c, double *loc_d)
{
    double tmp = 0;
    int i = 0;

    loc_c[0] = c[0] / b[0];  // c[0] - всегда 1  (-1)есть
    loc_d[0] = d[0] / b[0]; // d[0] - всегда 0 !!!

    for (i = 1; i<n; i++){
        tmp = b[i] - loc_c[i-1]*a[i];
        loc_c[i] = c[i] / tmp;
        loc_d[i] = (d[i*step] - loc_d[i-1]*a[i]) / tmp;
    }

    x[n*step-step] = loc_d[n - 1];
    for (i = (n-2); i>=0; i--)
        x[i*step] = loc_d[i] - loc_c[i]*x[i*step + step];

    return;
}

} // namespace Bpde
