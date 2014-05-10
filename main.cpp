#include <iostream>

//#include "defs.h"
#include "BArea.h"
#include "bsolveromp.h"


#include <cmath>
// temporary main file of srw project

int main(int argc, char **argv) { 

    using namespace std;
  
    int I = 100, J = 100, T = 2;
    int n = (I- 2)*(J -2);

    using namespace Boussinesq;

    BArea* area = new BArea("xyz0.txt");
    //комментарий на русском. тест
//    BArea* area           = new BArea(1, 1, 1, 300, 300, 100000);
    BSolverOmp* solver = new BSolverOmp(area, 10e-5, 5);
    
    solver->solve();

    return 0;
}
