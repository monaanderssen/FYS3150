#include "system.h"
#include "velocityverlet.h"
#include "lennardjones.h"
#include "statisticssampler.h"
#include "unitconverter.h"
#include "math/random.h"
#include <iostream>
#include "armadillo"

using namespace std;

System::System()
{

}

System::~System()
{
    for(Atom *atom : m_atoms) {
        delete atom;
    }
    m_atoms.clear();
}

void System::applyPeriodicBoundaryConditions() {
    for(Atom *atom : m_atoms) {
        for(int i=0; i<3; i++){
            if (atom->position[i] <  -m_systemSize[i] * 0.5) atom->position[i] = atom->position[i] + m_systemSize[i];
            if (atom->position[i] >=  m_systemSize[i] * 0.5) atom->position[i] = atom->position[i] - m_systemSize[i];
        }
    }
}

void System::removeTotalMomentum() {
    // Find the total momentum and remove momentum equally on each atom so the total momentum becomes zero.
    statisticsSampler.sampleTotalMomentum(*this);
    vec3 m_totalMomentum = statisticsSampler.totalMomentum();
    //m_totalMomentum = statisticsSampler.sampleTotalMomentum(*this);
    cout << m_totalMomentum << endl;
    for(Atom *atom : m_atoms){
        for(int i=0; i<3; i++){
            //cout << atom -> velocity[i];
            atom->velocity[i] -= m_totalMomentum[i]/m_numberOfAtoms/atom->mass();
            //cout << atom -> velocity[i] << endl;
        }
    }
    statisticsSampler.sampleTotalMomentum(*this);
    m_totalMomentum = statisticsSampler.totalMomentum();
    cout << m_totalMomentum << endl;
}

void System::createFCCLattice(int numberOfUnitCellsEachDimension, double latticeConstant, double temperature) {
    // You should implement this function properly. Right now, 100 atoms are created uniformly placed in the system of size (10, 10, 10).

    for(int i=0; i<100; i++) {
        Atom *atom = new Atom(UnitConverter::massFromSI(6.63352088e-26));
        double x = Random::nextDouble(0, 10); // random number in the interval [0,10]
        double y = Random::nextDouble(0, 10);
        double z = Random::nextDouble(0, 10);
        atom->position.set(x,y,z);
        atom->resetVelocityMaxwellian(temperature);
        m_atoms.push_back(atom);
    }
    cout << m_atoms.size() << endl;
    m_numberOfAtoms = m_atoms.size();
    setSystemSize(vec3(10, 10, 10)); // Remember to set the correct system size!
}

void System::createFCCLatticeCrystalStructure(int numberOfUnitCellsEachDimension, double latticeConstant, double temperature) {
    double b = 1.0;
    arma::vec x = {0, b/2.0, 0, b/2};
    arma::vec y = {0., b/2., b/2., 0.};
    arma::vec z = {0.,0.,b/2.,b/2.};
    for(int i=0; i<x.size(); i++) {
        Atom *atom = new Atom(UnitConverter::massFromSI(6.63352088e-26));

        atom->position.set(x[i],y[i],z[i]);
        atom->resetVelocityMaxwellian(temperature);
        m_atoms.push_back(atom);
    }
    cout << m_atoms.size() << endl;
    m_numberOfAtoms = m_atoms.size();
    setSystemSize(vec3(10, 10, 10)); // Remember to set the correct system size!
}

void System::calculateForces() {
    for(Atom *atom : m_atoms) {
        atom->resetForce();
    }
    m_potential.calculateForces(*this); // this is a pointer, *this is a reference to this object
}

void System::step(double dt) {
    m_integrator.integrate(*this, dt);
    m_steps++;
    m_time += dt;
}
