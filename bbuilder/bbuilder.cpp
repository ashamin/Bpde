#include "bbuilder.h"

namespace Bpde
{

BSolverBuilder* BSolverBuilder::instance = NULL;

BSolverBuilder::BSolverBuilder()
    :area(NULL)
{
}

BSolverBuilder* BSolverBuilder::getInstance()
{
    if (instance == NULL)
        return new BSolverBuilder();
    return instance;
}

BSolver* BSolverBuilder::getSolver(std::string file,
        ParallelizationMethod::ParallelizationMethod pMethod, int threadsNum)
{
    if (area != NULL)
        delete area;
    area = new BArea(file);
    if (pMethod == ParallelizationMethod::OPENMP)
        return new BSolverOmp(*area, threadsNum);
    return new BSolverOmp(BArea(file), 1);
}

}
