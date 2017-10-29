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

using namespace std;

int main()
{
    int IntegrationPoints;  // No. of integration points
    double FinalTime;       // End time of calculation
    int Dimension;           // No. of spatial dimensions

    //Removing files
    remove("/Users/monaanderssen/Documents/FYS3150/ProjectsInFYS3150/PlanetPositionEulerResults.txt");
    remove("/Users/monaanderssen/Documents/FYS3150/ProjectsInFYS3150/PlanetEnergiesEulerResults.txt");
    remove("/Users/monaanderssen/Documents/FYS3150/ProjectsInFYS3150/PlanetPositionVerletResults.txt");
    remove("/Users/monaanderssen/Documents/FYS3150/ProjectsInFYS3150/PlanetEnergiesVerletResults.txt");
    remove("/Users/monaanderssen/Documents/FYS3150/ProjectsInFYS3150/MercuryPostionResults.txt");

        cout << "Earth-Sun binary system" << endl;
        Dimension = 2;

        FinalTime =250.;
        IntegrationPoints = 1e5*FinalTime;
        //IntegrationPoints = 7*100*3600*360; // # integrationpoints to find the perihelion

        double epsilon = 0.0;

        //double TimeStep = FinalTime/IntegrationPoints;
        double x[3],v[3];  // positions and velocities

//        planet Sun(1.,0.,0.,0.,0.,0.,0.); // The sun without any movement
        planet Sun(1.,2.187003065211543E-03,5.768166559108312E-03,-1.294147734354897E-04,-0.00192864,0.00199457, 0.0000454587); // The sun with real values from NASA
//        planet Sun(1.,0.,0.,0.,0.,-2.636e-3,0.); // The sun with velocity

//        planet Earth(0.000003,1.,0.0,0.0,0.0,2*M_PI,0.); //Earth's values from table in the text
//        planet Earth(0.000003,1.,0.0,0.0,0.0,8.9,0.); //Earth, testing escape velocity
        planet Earth(0.000003,8.589987108796383E-01,5.110680605545673E-01,-1.568623415833688E-04,-3.29862,5.39202,-0.0000709728); //Earth with values from NASA

//        planet Jupiter(0.95e-3, 5.2, 0.0, 0.0, 0.0, 2*M_PI/(sqrt(5.2)), 0.0); //Jupiter's values from table in the text
        planet Jupiter(0.95e-3, -4.556745348155565E+00, -2.963008457339381E+00, 1.142108603087190E-01, 1.47024, -2.17987, -0.0238383); // Values from NASA
        planet Mars(3.3e-7, -1.590488403238053E+00, 4.879788693373922E-01, 4.906264799695509E-02, -1.29166, -4.45322, -0.0616634); // Values from NASA
        planet Venus(2.45e-6, -6.915411411024813E-01, 1.907707656342085E-01, 4.244080591200550E-02, -1.94245, -7.16966, 0.0136515); // Values from NASA
        planet Saturn(2.75e-4, -3.158467085324504E-01, -1.005065028034512E+01, 1.873222298678130E-01, 1.92497, -0.070396, -0.0754838); // Values from NASA

//        planet Mercury(1.65e-7, 0.3075, 0., 0., 0.,12.44, 0. ); // Values for Mercury to find the perihelion
        planet Mercury(1.65e-7, -2.139370590565288e-1, -4.028814669327753e-1, -1.369419923866817e-2, 7.01136, -4.30608, -0.995377); // Values from NASA

        planet Uranus(4.4e-5, 1.784724616993690E+01, 8.833225342557650E+00,-1.984072076144117E-01, -0.647764, 1.22053, 0.0129038);// Values from NASA
        planet Neptun(0.515e-4, 2.862016046630078E+01, -8.800209679340895E+00, -4.783572794559496E-01, 0.329248, 1.10259, -0.0303563); // Values from NASA
        planet Pluto(0.655e-8, 10.56835116967967, -3.171023274407168E+01, 3.361865390902906E-01, 1.10908, 0.121938, -0.330858); // Values from NASA

        solver binary_vv(5.0);

        binary_vv.setFileWriting(true,1000);

        binary_vv.add(Sun);
        binary_vv.add(Earth);
        binary_vv.add(Jupiter);
        binary_vv.add(Mars);
        binary_vv.add(Venus);
        binary_vv.add(Saturn);
        binary_vv.add(Mercury);
        binary_vv.add(Uranus);
        binary_vv.add(Neptun);
        binary_vv.add(Pluto);

        // Calling the Euler function
//        binary_vv.Euler(Dimension, IntegrationPoints, FinalTime, epsilon);

        // Calling the general Verlet function
        binary_vv.VelocityVerlet(Dimension,IntegrationPoints,FinalTime,epsilon);

        // Calling the Verlet function for Mercury
//        binary_vv.VelocityVerletMercury(Dimension,IntegrationPoints,FinalTime,epsilon, Sun, Mercury);

        for(int j = 0; j < Dimension;j++){
            x[j] = binary_vv.all_planets[0].position[j];
            v[j] = binary_vv.all_planets[0].velocity[j];
        }

    return 0;
}
