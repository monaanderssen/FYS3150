#include "system.h"
#include "statisticssampler.h"
#include "lennardjones.h"
#include <iostream>
#include "math.h"
#include <iomanip>
#include "unitconverter.h"

using std::ofstream; using std::cout; using std::endl;
using namespace std;


StatisticsSampler::StatisticsSampler()
{

}

void StatisticsSampler::saveToFile(System &system)
{
    // Save the statistical properties for each timestep for plotting etc.
    // First, open the file if it's not open already
    if(!m_file.good()) {
        m_file.open("/Users/monaanderssen/Documents/FYS3150/FYS3150/Project5/results/statistics.txt", ofstream::out);
        // If it's still not open, something bad happened...
        if(!m_file.good()) {
            cout << "Error, could not open statistics.txt" << endl;
            exit(1);
        }
    }
    //cout << system.steps() << endl;
    m_file << setw(20) << system.steps() <<
            setw(20) << system.time() <<
            setw(20) << UnitConverter::temperatureToSI(temperature() ) <<
            setw(20) << m_kineticEnergy <<
            setw(20) << m_potentialEnergy <<
            setw(20) << totalEnergy() <<
            setw(20) << m_diffusionConstant << endl;
}

void StatisticsSampler::sample(System &system)
{
    // Here you should measure different kinds of statistical properties and save it to a file.
    sampleKineticEnergy(system);
    samplePotentialEnergy(system);
    sampleTemperature(system);
    sampleDiffusionConstant(system);
    //sampleDensity(system);
    saveToFile(system);
}

void StatisticsSampler::sampleKineticEnergy(System &system)
{
    m_kineticEnergy = 0; // Remember to reset the value from the previous timestep
    for(Atom *atom : system.atoms()) {
        m_kineticEnergy += 0.5*atom->mass()*atom->velocity.lengthSquared();
    }
}

void StatisticsSampler::samplePotentialEnergy(System &system)
{
    m_potentialEnergy = system.potential().potentialEnergy();
}

void StatisticsSampler::sampleTotalMomentum(System &system){
    m_totalMomentum.zeros();
    for(Atom *atom : system.atoms()) {
        for(int i=0; i<3; i++){
            m_totalMomentum[i] += atom->velocity[i]*atom->mass();
        }
    }
    std::cout << m_totalMomentum << std::endl;
}

void StatisticsSampler::sampleTemperature(System &system)
{
    m_temperature = (2.0/3.0*m_kineticEnergy/system.atoms().size());
}

void StatisticsSampler::sampleDensity(System &system, int N_x, int N_y, int N_z){
    m_density = 0.0;
    double m_b = system.b();
    for(Atom *atom : system.atoms()) {
        m_density += (atom->mass())/(m_b*m_b*m_b*N_x*N_y*N_z);
    }
}

void StatisticsSampler::sampleDiffusionConstant(System &system){
    double temp_diffusion = 0;
    for(int i=0; i < system.numberOfAtoms(); i++){
        Atom* atom = system.atoms()[i];
        temp_diffusion += (atom->m_distanceBeforePBC * atom->m_distanceBeforePBC)/(6*system.time());
    m_diffusionConstant = temp_diffusion/system.numberOfAtoms();
    }
}
