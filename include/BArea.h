#ifndef B_AREA
#define B_AREA

#include <string>

namespace Bpde
{

class BArea {
public:
    /**
     * @brief BArea
     * @param sizeX размер области по x
     * @param sizeY размер области по y
     * @param destTime
     * @param I количество узлов по x
     * @param J количество узлов по y
     * @param T количество разбиений по времени. определяет шаг по времени dt
     */
    BArea(double sizeX, double sizeY, double destTime, int I, int J, double T);
    BArea(std::string file);
    BArea(const BArea& area);
    ~BArea();

    double answer(double x, double y, double t);

    virtual double V(double x, double y, double t);

    double* H;
    double* x;
    double* y;
    int I, J, T;

    double zc, zf;//= 160, zf = 90;
    double mu1, mu2;//= 0.16, mu2 = 0.16;
    double kx, ky;//= (double)1/(3600*24), ky = (double)1/(3600*24);

//private:
    double dt;
};

} // namespace Bpde

#endif
