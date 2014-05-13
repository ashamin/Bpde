#ifndef BSOLVER_OMP
#define BSOLVER_OMP

#include "BArea.h"
#include "bsolver.h"
#include "omp.h"

#include <iostream>
#include "math.h"

namespace Bpde{

/// TODO: move theese variables to __bpde_omp namespace
namespace __bpde_omp {
    static double *H, *Ha;
    static double *x, *y;
    static double *V;
    static double *mu;
    static double *dx_d, *dx_l, *dx_u, *dy_d, *dy_l, *dy_u;
    static double *loc_c, *loc_d, *b;
    static int I, J, T;
    static double dt;
}
///////////


// z ceiling and z floor
const double zc = 160, zf = 90;
const double mu1 = 0.16, mu2 = 0.16;
const double kx = (double)1/(3600*24), ky = (double)1/(3600*24);

inline void getMu(double* mu, double H)
{
    *mu = (H >= zc)?mu1:mu2;
}

inline double Tx(double H)
{
    return (H >= zc)?kx*(zc - zf):((H < zf)?0:kx*(H - zf));
}

inline double Ty(double H)
{
    return (H >= zc)?ky*(zc - zf):((H < zf)?0:ky*(H - zf));
}

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

class BSolverOmp : public BSolver
{
public:
    BSolverOmp(BArea* area, double epsilon, int maxit);
    virtual ~BSolverOmp();

	double* solve();

	// H computation area
	double* H;
	// approximation of H
	double* Ha;

	double* b;
	
	double epsilon;
	int maxit;

	int n;
	double* V;
	
	double exec_time();
	int it_num();

	BArea* area;

    double* x;
    double* y;

	double t;
    double dt;

private:
    inline void prepareIteration();

	// diags of X differential operator
	double* dx_d;
	double* dx_l;
	double* dx_u;

	// diags of Y differentioal operator
	double* dy_d;
	double* dy_l;
	double* dy_u;

	double* loc_c;
	double* loc_d;

	//coefficient near time differential
	double* mu;

	int I, J;

    int iterations;
};

}

#endif
