#ifndef LSOLVER_H
#define LSOLVER_H


/**
 * Base class for solvers of sparse linear systems.
 * \author ashamin
 */
class BSolver{
public:
  double* solve();
  double exec_time();
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


#endif
