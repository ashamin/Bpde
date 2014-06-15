#include "BLogger/BLoggerState.h"

namespace Bpde
{

BLoggerStateEnabled::BLoggerStateEnabled(std::basic_ostream<char> *stream)
{
    log = stream;
}

void BLoggerStateEnabled::logMatrix(std::string name, double *var, int xSz, int ySz)
{
    *log << "RESULT__" << name << ':' << std::endl;

    for (int j = 0; j<ySz; j++){
        for (int i = 0; i<xSz; i++)
            *log << var[j*xSz + i] << ' ';
            *log << std::endl;
    }
    *log << std::endl;
}

void BLoggerStateEnabled::logVector(std::string name, double *var, int size)
{
    *log << "RESULT__" << name << ':' << std::endl;
    for (int i = 0; i<size; i++)
        *log << var[i] << ' ';
    *log << std::endl;
}

void BLoggerStateEnabled::log3DiagMatrix(std::string name, double *lower, double *main,
        double *upper, int size, int sz)
{
    *log << "RESULT__" << name << ':' << std::endl;
    int k = 0;
    double eps = .00005;
   // int sz = (int)sqrt(size);
    for (int i = 0; i<size; i++){
        for (int j = 0; j<size; j++){
            if (j == k-1)
                *log << ((fabs(lower[k])>eps)?lower[k]:0);

            else if (j == k)
                *log << ((fabs(main[k])>eps)?main[k]:0);

            else if (j == k+1)
                *log << ((fabs(upper[k])>eps)?upper[k]:0);
            else
                *log << 0;
            *log << " ";
            //vertical divider
            if ((j+1)%sz==0) *log << "| ";
        }
        k++;
        *log << std::endl;
        // horizontal divider
        if ((i+1)%sz==0){
            for (int p = 0; p<size; p++)
                *log << "___";
            *log << std::endl;
        }
    }
}

} // namespace Bpde
