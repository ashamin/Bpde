#include "gtest/gtest.h"

#include "BBuilder.h"

TEST (QualityTests, test1)
{
//    SCOPED_TRACE("test1");
//    double zf = 90, zc =160;
//    using namespace Bpde;
//    BSolver* solver = BSolverBuilder::getInstance()->
//            getSolver("xyz0.txt", ParallelizationMethod::OPENMP);
//    solver->solve();
//    int I = __bpde_omp::I;
//    int J = __bpde_omp::J;
//    for (int i = 1; i<I-1; i++){
//        for (int j = 1; j<J-1; j++){
//            if (__bpde_omp::H[j*I+i] < zf || __bpde_omp::H[j*I+i] > zc)
//                FAIL() << "Not all nodes of grid has values between floor and ceining"
//                       << "first entry with coordinates (" << __bpde_omp::x[i] << ";"
//                       << __bpde_omp::y[j];
//            if (__bpde_omp::Ha[j*I+i] > 10e-9)
//                FAIL() << "Derivative did not converged to values of 10e-9 order";
//        }
//    }
//    delete solver;
}
