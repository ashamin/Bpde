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
//            getSolver("xyz0.txt", ParallelizationMethod::OPENMP, 1);
//    solver->solve();
//    std::cout << solver->exec_time() << std::endl;

//    delete solver;

//    solver = BSolverBuilder::getInstance()->
//                getSolver("xyz0.txt", ParallelizationMethod::OPENMP, 2);
//    solver->solve();
//    std::cout << solver->exec_time() << std::endl;

//    delete solver;



    std::vector<cl::Platform> platforms;
    std::vector<cl::Device> devices;

    try{
        cl::Platform::get(&platforms);
        platforms[0].getDevices(CL_DEVICE_TYPE_CPU, &devices);
    }
    catch(...) {
    }


    BArea area("xyz0.txt");
    BSolver* solver = new BSolverOcl(area, devices[0], 1);
    solver->solve();
    std::cout << solver->exec_time() << std::endl;


    return 0;
}
