#ifndef BBUILDER_H
#define BBUILDER_H

#include "BSolverOmp.h"

namespace Bpde
{

class BSolverBuilder
{
public:
    static BSolverBuilder* getInstance();

    BSolver* getSolver(std::string file, ParallelizationMethod::ParallelizationMethod pMethod =
            ParallelizationMethod::OPENMP, int threadsNum = omp_get_max_threads());

private:
    BSolverBuilder();
    static BSolverBuilder* instance;

    BArea* area;
};

} // namespace Bpde

#endif // BBUILDER_H