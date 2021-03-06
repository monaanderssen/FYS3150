#include "math/random.h"
#include "lennardjones.h"
#include "velocityverlet.h"
#include "system.h"
#include "statisticssampler.h"
#include "atom.h"
#include "io.h"
#include "unitconverter.h"
#include <iostream>
#include <iomanip>

using namespace std;


int main(int numberOfArguments, char **argumentList)
{
    remove("/Users/monaanderssen/Documents/FYS3150/Project5/results/diffusion_different_temperatures.txt");
    //remove("/home/pederbh/UiO/FYS4150/FYS3150/Project5/results/diffusion_different_temperatures.txt");
    int numberOfUnitCells = 5;
    int N_x = 4; int N_y= 4; int N_z = 4;
    int max_time;
    double input_temperature, max_temperature;
    if(numberOfArguments > 1) max_time = atoi(argumentList[1]);
    if(numberOfArguments > 2) input_temperature = atof(argumentList[2]);
    if(numberOfArguments > 3) max_temperature = atof(argumentList[3]);

    double latticeConstant = UnitConverter::lengthFromAngstroms(5.26); // measured in angstroms

    // If a first argument is provided, it is the number of unit cells
    //if(numberOfArguments > 1) numberOfUnitCells = atoi(argumentList[1]);
    //if(numberOfArguments > 1) max_time = atoi(argumentList[1]);
    // If a second argument is provided, it is the initial temperature (measured in kelvin)
    //if(numberOfArguments > 2) initialTemperature = UnitConverter::temperatureFromSI(atof(argumentList[2]));
    //if(numberOfArguments > 2) input_temperature = atof(argv[2]);
    // If a third argument is provided, it is the lattice constant determining the density (measured in angstroms)
    //if(numberOfArguments > 3) latticeConstant = UnitConverter::lengthFromAngstroms(atof(argumentList[3]));

    double dt = UnitConverter::timeFromSI(1e-15); // Measured in seconds.

    cout << "One unit of length is " << UnitConverter::lengthToSI(1.0) << " meters" << endl;
    cout << "One unit of velocity is " << UnitConverter::velocityToSI(1.0) << " meters/second" << endl;
    cout << "One unit of time is " << UnitConverter::timeToSI(1.0) << " seconds" << endl;
    cout << "One unit of mass is " << UnitConverter::massToSI(1.0) << " kg" << endl;
    cout << "One unit of temperature is " << UnitConverter::temperatureToSI(1.0) << " K" << endl;

    for(float current_temperature = input_temperature; current_temperature < max_temperature; current_temperature = current_temperature + 10){
        double initialTemperature = UnitConverter::temperatureFromSI(current_temperature); // measured in Kelvin
        System system;
        system.createFCCLatticeCrystalStructure(numberOfUnitCells, latticeConstant, initialTemperature, N_x, N_y, N_z);
        system.potential().setEpsilon(UnitConverter::temperatureFromSI(119.8)); //119.8 or 1?
        system.potential().setSigma(UnitConverter::lengthFromAngstroms(3.405));

        system.removeTotalMomentum();
        double diffusion_average = 0;
        double temperature_average = 0;
        double average_counter = 0;

        StatisticsSampler statisticsSampler;
        statisticsSampler.sampleTemperature(system);
        statisticsSampler.m_initialTemperature = initialTemperature;

        //system.potential().setEpsilon(UnitConverter::temperatureFromSI(119.8));
        statisticsSampler.sampleDensity(system, N_x, N_y, N_z);
        cout << "Density of the system: " << statisticsSampler.density() << endl;

        IO movie("/Users/monaanderssen/Documents/FYS3150/Project5/results/movie.xyz");
        //IO movie("/home/pederbh/UiO/FYS4150/FYS3150/Project5/results/movie.xyz"); // To write the state to file
        cout << "max" << max_time << endl;
        cout << setw(20) << "Timestep" <<
                setw(20) << "Time" <<
                setw(20) << "Temperature" <<
                setw(20) << "KineticEnergy" <<
                setw(20) << "PotentialEnergy" <<
                setw(20) << "TotalEnergy" <<
                setw(20) << "Diffusion constant" << endl;
        for(int timestep=0; timestep<(int)max_time; timestep++) {
            system.step(dt);
            statisticsSampler.sample(system);
            if( timestep % 1000 == 0 ) {
                // Print the timestep every 100 timesteps
                cout << setw(20) << system.steps() <<
                        setw(20) << system.time() <<
                        setw(20) << UnitConverter::temperatureToSI(statisticsSampler.temperature() ) <<
                        setw(20) << statisticsSampler.kineticEnergy() <<
                        setw(20) << statisticsSampler.potentialEnergy() <<
                        setw(20) << statisticsSampler.totalEnergy() <<
                        setw(20) << statisticsSampler.diffusionConstant() << endl;
            }
            //if( timestep % 10 == 0 )  movie.saveState(system);
            movie.saveState(system);

            if( timestep > (max_time - 50) ){ //last 100 iterations, write diffusion constant to file
                diffusion_average += statisticsSampler.diffusionConstant();
                temperature_average += statisticsSampler.temperature();
                average_counter += 1;
                if( timestep == (max_time - 1)){
                    string path= string("/Users/monaanderssen/Documents/FYS3150/Project5/results/diffusion_different_temperatures.txt");
                    //string path= string("/home/pederbh/UiO/FYS4150/FYS3150/Project5/results/diffusion_different_temperatures.txt");
                    diffusion_average /= average_counter;
                    temperature_average /= average_counter;
                    ofstream diffusionFile;
                    diffusionFile.open(path, std::ios::app);
                    diffusionFile << setw(20) << UnitConverter::temperatureToSI(temperature_average) << " " << UnitConverter::diffusionToSI(diffusion_average) << endl;
                    diffusionFile.close();
                    cout << "Wrote to diffusion file." << endl;
                }
            }
        }

        movie.close();
    }
    return 0;
}
