#include <iostream>

//#include "defs.h"
#include "BArea.h"
#include "BSolverOmp.h"
#include "BSolverOcl.h"
#include "BBuilder.h"

#include "CL/cl.h"

#include <cmath>

#include <iostream>
#include <stdio.h>

// temporary main file of srw project
using namespace Bpde;
int main(int argc, char **argv) {

    using namespace std;

    int I = 100, J = 100, T = 2;
    int n = (I- 2)*(J -2);

    using namespace Bpde;

// usage 1
//    BArea area("xyz0.txt");
//    BSolverOmp solver(area, omp_get_max_threads());
//    solver.solve();
//    std::cout << solver.exec_time() << std::endl;

// usage 2
//    BArea area("xyz0.txt");
//    BSolver* solver = new BSolverOmp(area, omp_get_max_threads());
//    solver->solve();
//    std::cout << solver->exec_time() << std::endl;

// maintain usage
//    BSolver* solver = BSolverBuilder::getInstance()->
//            getSolver("xyz0.txt", ParallelizationMethod::OPENMP);
//    solver->solve();
//    std::cout << solver->exec_time() << std::endl;

//    delete solver;

//    BSolver* solver = BSolverBuilder::getInstance()->
//                getSolver("xyz0.txt", ParallelizationMethod::OPENMP, 2);
//    solver->solve();
//    std::cout << solver->exec_time() << std::endl;

//    delete solver;

// opencl usage 1
//    std::vector<cl::Platform> platforms;
//    std::vector<cl::Device> devices;

//    try{
//        cl::Platform::get(&platforms);
//        platforms[0].getDevices(CL_DEVICE_TYPE_CPU, &devices);
//    }
//    catch(...) {
//    }


//    BArea area("xyz0.txt");
//    BSolver* solver = new BSolverOcl(area, devices);
//    solver->solve();
//    std::cout << solver->exec_time() << std::endl;

// maintaint opencl usage
//    BSolver* solver = BSolverBuilder::getInstance()->
//            getSolver("xyz0.txt", ParallelizationMethod::OPENCL);
//    solver->solve();

//    std::cout << solver->exec_time() << std::endl;
//    delete solver;


//    for (int i = 0; i<10; i++)
//    {

        BSolver* solver = BSolverBuilder::getInstance()->
                getSolver("xyz0.txt", ParallelizationMethod::OPENMP, 2);
        solver->solve();
        std::cout << solver->exec_time() << std::endl;

        delete solver;
//    }


    return 0;
}
