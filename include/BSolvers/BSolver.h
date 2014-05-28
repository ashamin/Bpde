#ifndef BSOLVER_H
#define BSOLVER_H

#include <iostream>
#include <cmath>

#include "BArea.h"

namespace Bpde {

namespace ParallelizationMethod
{
enum ParallelizationMethod
{
    NONE = 0,
    OPENMP = 1,
    OPENCL = 2
};
} // namespace ParallelizationMethod

// z ceiling and z floor

inline void log_matrix(char *name, double* var, int xSz, int ySz)
{
    std::cout << "RESULT__" << name << ':' << std::endl;

    for (int j = 0; j<ySz; j++){
        for (int i = 0; i<xSz; i++)
            std::cout << var[j*xSz + i] << ' ';
            std::cout << std::endl;
    }
    std::cout << std::endl;
}

inline void log_vector(char *name, double* var, int size)
{
    std::cout << "RESULT__" << name << ':' << std::endl;
    for (int i = 0; i<size; i++)
        std::cout << var[i] << ' ';
    std::cout << std::endl;
}

inline void log_diags_as_3dmatrix(char *name,
        double* lower,
        double* main,
        double* upper, int size, int sz)
{
    std::cout << "RESULT__" << name << ':' << std::endl;
    int k = 0;
    double eps = .00005;
   // int sz = (int)sqrt(size);
    for (int i = 0; i<size; i++){
        for (int j = 0; j<size; j++){
            if (j == k-1)
                std::cout << ((fabs(lower[k])>eps)?lower[k]:0);

            else if (j == k)
                std::cout << ((fabs(main[k])>eps)?main[k]:0);

            else if (j == k+1)
                std::cout << ((fabs(upper[k])>eps)?upper[k]:0);
            else
                std::cout << 0;
            std::cout << " ";
            //vertical divider
            if ((j+1)%sz==0) std::cout << "| ";
        }
        k++;
        std::cout << std::endl;
        // horizontal divider
        if ((i+1)%sz==0){
            for (int p = 0; p<size; p++)
                std::cout << "___";
            std::cout << std::endl;
        }
    }
}

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
    virtual void addExtraIterations(int its) = 0;
    virtual void setTimeStep(double dt) = 0;
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
