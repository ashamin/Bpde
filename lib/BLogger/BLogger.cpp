#include "BLogger/BLogger.h"

namespace Bpde
{

BLogger* BLogger::instance = NULL;

void BLogger::enable()
{
    delete state;
    state = new BLoggerStateEnabled(stream);
}

void BLogger::disable()
{
    delete state;
    state = new BLoggerStateDisabled();
}

void BLogger::setStream(std::basic_ostream<char> *stream)
{
    this->stream = stream;
    if (dynamic_cast<BLoggerStateEnabled*>(state)) {
        delete state;
        state = new BLoggerStateEnabled(stream);
    }
}

void BLogger::logMatrix(std::string name, double *var, int xSz, int ySz)
{
    state->logMatrix(name, var, xSz, ySz);
}

void BLogger::log3DiagMatrix(std::string name, double *lower, double *main, double *upper,
        int size, int sz)
{
    state->log3DiagMatrix(name, lower, main, upper, size, sz);
}

void BLogger::logVector(std::string name, double *var, int size)
{
    state->logVector(name, var, size);
}

BLogger::BLogger()
    : state(NULL), stream(&std::cout)
{
    state = new BLoggerStateDisabled();
}

} // namespace Bpde
