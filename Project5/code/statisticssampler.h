#ifndef STATISTICSSAMPLER_H
#define STATISTICSSAMPLER_H
#include <fstream>
#include "math/vec3.h"

class System; // Promise the compiler that this is a class even though we haven't included system.h here

class StatisticsSampler
{
private:
    std::ofstream m_file;
    double m_kineticEnergy = 0;
    double m_potentialEnergy = 0;
    vec3 m_totalMomentum = {0,0,0};
    double m_temperature = 0;
    double m_density = 0;
    double m_diffusionConstant;
    double m_r2;
    double m_TRatio;
public:
    double m_initialTemperature;
    StatisticsSampler();
    void saveToFile(System &system);
    void sampleTRatio(System &system);
    void saveDiffusionDifferentTemperatures(System &system);
    void sample(System &system);
    void sampleKineticEnergy(System &system);
    void samplePotentialEnergy(System &system);
    void sampleTotalMomentum(System &system);
    void sampleTemperature(System &system);
    void sampleDensity(System &system, int N_x, int N_y, int N_z);
    void sampleDiffusionConstant(System &system);
    double kineticEnergy() { return m_kineticEnergy; }
    double potentialEnergy() { return m_potentialEnergy; }
    double totalEnergy() { return m_kineticEnergy+m_potentialEnergy; }
    double temperature() { return m_temperature; }
    double density() { return m_density; }
    double diffusionConstant() { return m_diffusionConstant; }
    vec3 totalMomentum(){ return m_totalMomentum; }
};
#endif
