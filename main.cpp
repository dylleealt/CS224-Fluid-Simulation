#include "Solver3D.h"

int main(){
    Solver3D *solver = new Solver3D();
    solver->init();
    solver->reset();

    // while (1){
    for (int i = 0; i < 100; ++i){
       solver->vStep();
       solver->sStep();
    }

    return 0;
}
