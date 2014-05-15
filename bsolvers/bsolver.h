#ifndef BSOLVER_H
#define BSOLVER_H

namespace Bpde {

namespace ParallelizationMethod
{
enum ParallelizationMethod{
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
    BSolver(){}
    virtual ~BSolver(){}
    virtual double* solve() = 0;
    virtual double exec_time() = 0;
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
