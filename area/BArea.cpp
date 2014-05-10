#include "BArea.h"
#include <cmath>

#include "stdlib.h"
#include <iostream>
#include <fstream>

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

    area >> I;
    area >> hx;
    area >> J;
    area >> hy;
    area >> T;
    area >> dt;
    double tmp;

//    test << I << " " << hx << " " << J << " " << hy << " " << T << " " << dt << std::endl;

    H = new double((I+2)*(J+2));

    for (int i = 1; i<I+1; i++){
        for (int j = 1; j<J+1; j++){
            area >> tmp;
            area >> tmp;
            area >> H[i+(I+2)*j];
//            test << (i-1)*hx << " " << (j-1)*hy << " " << H[j+(J+2)*i] << std::endl;
        }
    }

    area.close();

    for (int j = 0; j<J+2; j++){
        for (int i = 0; i<I+2; i++)
           test <<  H[i+(I+2)*j] << "\t";
        test << std::endl;
    }
    test.close();
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
