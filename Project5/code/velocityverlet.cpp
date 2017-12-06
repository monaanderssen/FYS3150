#include "velocityverlet.h"
#include "system.h"
#include "atom.h"
#include "math.h"

void VelocityVerlet::integrate(System &system, double dt)
{
    if(m_firstStep) {
        system.calculateForces();
        m_firstStep = false;
    }

    for(Atom *atom : system.atoms()) {
        atom->velocity += atom->force*0.5*dt/atom->mass();
        atom->position += atom->velocity*dt;
        double dx = atom->position[0] - atom->initialPosition[0];
        double dy = atom->position[1] - atom->initialPosition[1];
        double dz = atom->position[2] - atom->initialPosition[2];
        double r = sqrt(dx*dx + dy*dy + dz*dz);
        atom->m_distanceBeforePBC = r;
    }

    system.applyPeriodicBoundaryConditions();
    system.calculateForces(); // New positions, recompute forces

    for(Atom *atom : system.atoms()) {
        atom->velocity += atom->force*0.5*dt/atom->mass();
    }
}
