#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <cstring>
#include <fstream>
#include <random>
#include <chrono>
#include <time.h>
#include "planet.h"
#include "solver.h"
#include <cstdio>
#include "armadillo"

using namespace std;

int main()
{

}

void EnergyAndMagnetization(){
    double E, M, E_2, M_2;

    for(int i=0;i<N;i++){
        E = (1/N)*(E(i));
        M = (1/N)*(M(i));
        E_2 = (1/N)*(E(i)*E(i));
        M_2 = (1/N)*(M(i)*M(i));
    }
}
