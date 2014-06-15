#ifndef BLOGGERSTATE_H
#define BLOGGERSTATE_H

#include <string>
#include <iostream>
#include <cmath>

namespace Bpde
{

class BLoggerState
{
public:
    BLoggerState(){}
    virtual ~BLoggerState(){}
    inline virtual void logMatrix(std::string name, double* var, int xSz, int ySz) = 0;
    inline virtual void log3DiagMatrix(std::string name, double* lower,
            double *main, double *upper, int size, int sz) = 0;
    inline virtual void logVector(std::string name, double* var, int size) = 0;
}; // BLoggerState

class BLoggerStateEnabled : public BLoggerState
{
public:
    BLoggerStateEnabled(){log = &std::cout;}
    BLoggerStateEnabled(std::basic_ostream<char>* stream);
    inline virtual ~BLoggerStateEnabled(){delete log;}
    inline virtual void logMatrix(std::string name, double* var, int xSz, int ySz);
    inline virtual void log3DiagMatrix(std::string name, double* lower,
            double *main, double *upper, int size, int sz);
    inline virtual void logVector(std::string name, double* var, int size);
private:
    std::basic_ostream<char>* log;
}; // BLoggerStateEnabled

class BLoggerStateDisabled : public BLoggerState
{
public:
    BLoggerStateDisabled(){}
    virtual ~BLoggerStateDisabled(){}
    inline virtual void logMatrix(std::string name, double* var, int xSz, int ySz){}
    inline virtual void log3DiagMatrix(std::string name, double* lower,
            double *main, double *upper, int size, int sz) {}
    inline virtual void logVector(std::string name, double* var, int size) {}
}; // BLoggerStateDisabled

} // namespace Bpde


#endif
