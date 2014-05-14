#include "BArea.h"
#include <cmath>

#include "stdlib.h"
#include <iostream>
#include <fstream>

namespace Bpde
{

BArea::BArea(double sizeX, double sizeY, double destTime, int I, int J, double T){
    this->I = I;
    this->J = J;
    this->T = T;
    hx = sizeX / (I - 1);
    hy = sizeY / (J - 1);
    dt = destTime / (T - 1);
}

BArea::BArea(std::string file)
{
    std::ifstream area;
    area.open(file.c_str());
    if(!area)
    {
        std::cerr << "Unable to open file " << file;
        exit(1);
    }

    std::ofstream test;
    test.open("olol.txt", std::ios::out);

    hx = 50;
    hy = 50;

    area >> I;
    area >> J;
    area >> T;
    area >> dt;

    I += 2;
    J += 2;

    H = new double[I*J];
    x = new double[I];
    y = new double[J];

    for (int i = 0; i<I*J; i++)
        H[i] = 0;
    x[0] = y[0] = x[I-1] = y[J-1] = 0.0;

    for (int i = 1; i<I-1; i++){
        for (int j = 1; j<J-1; j++){
            area >> x[i];
            area >> y[j];
            area >> H[i+I*j];
        }
    }

    area.close();

//    test << I << " " << J << " " << T << " " << dt << std::endl;

//    for (int i = 1; i<I+1; i++)
//        for (int j = 1; j<J+1; j++)
//            test << x[i] << " " << y[j] << " " << H[i+(I+2)*j] << std::endl;
//    for (int j = 0; j<J+2; j++){
//        for (int i = 0; i<I+2; i++)
//           test <<  H[i+(I+2)*j] << "\t";
//        test << std::endl;
//    }
//    test.close();
    for (int j = 0; j<J-1; j++){
            for (int i = 0; i<I-1; i++)
               test <<  H[i+I*j] << "\t";
            test << std::endl;
        }
        test.close();
}

BArea::BArea(const BArea &area)
{
    this->I = area.I;
    this->J = area.J;
    this->T = area.T;
    this->dt = area.dt;

    H = new double[I*J];
    x = new double[I];
    y = new double[J];

    hx = 50;
    hy = 50;

    for (int i = 0; i<I; i++)
        x[i] = area.x[i];
    for (int j = 0; j<J; j++)
        y[j] = area.y[j];
    for (int i = 0; i<I*J; i++)
        H[i] = area.H[i];
}

BArea::~BArea()
{
    delete[] H;
    delete[] x;
    delete[] y;
}

double BArea::answer(double x, double y, double t){
    return H[static_cast<int>(x/hx)+(I+2)*static_cast<int>(y/hy)];
//    return cos(x*3.14159265359)*cos(y*3.14159265359);
//    return cos(3.14159265359 * x);
//    return 1 + 3*y*y*y;
}

double BArea::V(double x, double y, double t){
    return 0;
//    return 2*cos(x*3.14159265359)*cos(y*3.14159265359) * (3.14159265359*3.14159265359);
//    return 3.14159265359 * 3.14159265359 * cos(3.14159265359 * x);
//    return -18*y;
}

} // namespace Bpde
