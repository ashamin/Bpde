Bpde
====

Bpde is demo software for solving Boussinesq partial differential equation. This equatuion is shown below Bpde gives interface for getting solvers with specific parameters which includes areas, dimentional and time steps, coefficients of Boussinesq equation, parallelization parameters. Also Bpde can be compiled to appication with GUI called BpdeGui. This GUI were desingned to demonstrate capabilities of bpde library and performance of parallel computations.

###Parallelization tuning
bpde library gives you interface to getting solver with specific characteristics. 
Bpde provides a factory class called BSolverBuilder. getSolver method takes 3 parameters and return configured solver of BSolver class:

1. Path to text file with area written in special format. 
2. Parallelization method
3. Optional parameter
  - Number of threads for OpenMP parallelization method
  - Vector of cl::Device (OpenCL class) for OpenCL parallelization method

Third parameter can be omitted. If so, it will have value of maximum avaliable number of CPU's in case of OpenMP method, and all OpenCL-devices of computational system in case of OpenCL.

```
Bpde::BSolver* solver = Bpde::BSolverBuilder::getInstance()->
    getSolver("xyz0.txt", ParallelizationMethod::OPENMP, 2);
```

Computations on all avaliable GPU devices.

```
std::vector<cl::Platform> platforms;
std::vector<cl::Device> tmp;
cl::Platform::get(&platforms);
for (std::vector<cl::Platform>::iterator it = platforms.begin();
    it!=platforms.end(); it++)
{
  (*it).getDevices(CL_DEVICE_TYPE_GPU, &tmp);
  std::copy(tmp.begin(), tmp.end(), std::back_inserter(devices));
}
Bpde::BSolver* solver = Bpde::BSolverBuilder::getInstance()->getSolver(
    sourceFileEdit->text().toStdString(), Bpde::ParallelizationMethod::OPENCL, devices)

```
