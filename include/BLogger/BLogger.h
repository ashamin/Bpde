#ifndef BLOGGER_H
#define BLOGGER_H

#include "BLoggerState.h"

class BLogger
{
public:
    BLogger* getInstance()
    {
        if (instance == NULL)
            instance = new BLogger();
        return instance;
    }
    void enable();
    void disable();
    void logMatrix(std::string name, double* var, int xSz, int ySz);
    void log3DiagMatrix(std::string name, double* lower, double *main, double *upper,
            int size, int sz);
    void logVector(std::string name, double* var, int size);
private:
    BLogger();
    static BLogger* instance;
    BLoggerState* state;

};

#endif // BLOGGER_H
