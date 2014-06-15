#ifndef BLOGGER_H
#define BLOGGER_H

#include "BLoggerState.h"
#include <iostream>

namespace Bpde
{

class BLogger
{
public:
    static BLogger* getInstance()
    {
        if (instance == NULL)
            instance = new BLogger();
        return instance;
    }
    void enable();
    void disable();
    void setStream(std::basic_ostream<char>* stream);
    void logMatrix(std::string name, double* var, int xSz, int ySz);
    void log3DiagMatrix(std::string name, double* lower, double *main, double *upper,
            int size, int sz);
    void logVector(std::string name, double* var, int size);
private:
    BLogger();
    virtual ~BLogger(){delete stream;}
    static BLogger* instance;
    BLoggerState* state;
    std::basic_ostream<char>* stream;
};

} // namespace Bpde

#endif // BLOGGER_H
