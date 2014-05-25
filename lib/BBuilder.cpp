#include "BBuilder.h"

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
    if (pMethod == ParallelizationMethod::OPENCL)
    {
        std::vector<cl::Platform> platforms;
        std::vector<cl::Device> devices;
        std::vector<cl::Device> tmp;

        try {
            cl::Platform::get(&platforms);
            for (std::vector<cl::Platform>::iterator it = platforms.begin();
                 it!=platforms.end(); it++)
            {
                (*it).getDevices(CL_DEVICE_TYPE_ALL, &tmp);
                std::copy(tmp.begin(), tmp.end(), std::back_inserter(devices));
            }
        }
        catch(...) {
        }
        return new BSolverOcl(*area, devices);
    }
    return new BSolverOmp(*area, 1);
}

BSolver* BSolverBuilder::getSolver(std::string file,
        ParallelizationMethod::ParallelizationMethod pMethod, std::vector<cl::Device>& devices)
{
    // only if opencl
    if (area != NULL)
        delete area;
    area = new BArea(file);
    if (pMethod == ParallelizationMethod::OPENCL)
        return new BSolverOcl(*area, devices);
}

} // namespace Bpde
