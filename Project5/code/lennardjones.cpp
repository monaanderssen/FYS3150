#include "lennardjones.h"
#include "system.h"
#include "math.h"

double LennardJones::potentialEnergy() const
{
    return m_potentialEnergy;
}

double LennardJones::sigma() const
{
    return m_sigma;
}

void LennardJones::setSigma(double sigma)
{
    m_sigma = sigma;
}

double LennardJones::epsilon() const
{
    return m_epsilon;
}

void LennardJones::setEpsilon(double epsilon)
{
    m_epsilon = epsilon;
}

void LennardJones::calculateForces(System &system)
{
    for(int j=0; j < system.numberOfAtoms(); j++){
        Atom* atomj = system.atoms()[j];
        atomj->resetForce();
    }

    double U, r_ij, V, Fx, Fy, Fz, dx, dy, dz, temp_exp, r2_ij, temp_exp2, m_eps24;
    m_potentialEnergy = 0; // Remember to compute this in the loop
    for(int j=0; j < system.numberOfAtoms(); j++){
        Atom* atomj = system.atoms()[j];
        for(int i=j+1; i < system.numberOfAtoms(); i++){
            Atom* atomi = system.atoms()[i];
            //Fx = Fy = Fz = 0;
            dx = atomj->position[0] - atomi->position[0];
            dy = atomj->position[1] - atomi->position[1];
            dz = atomj->position[2] - atomi->position[2];

            //Boundary conditions
            if (dx >   system.systemSize()[0] * 0.5) dx = dx - system.systemSize()[0];
            if (dx <= -system.systemSize()[0] * 0.5) dx = dx + system.systemSize()[0];


            if (dy >   system.systemSize()[1] * 0.5) dy = dy - system.systemSize()[1];
            if (dy <= -system.systemSize()[1] * 0.5) dy = dy + system.systemSize()[1];


            if (dz >   system.systemSize()[2] * 0.5) dz = dz - system.systemSize()[2];
            if (dz <= -system.systemSize()[2] * 0.5) dz = dz + system.systemSize()[2];

            r2_ij = dx*dx + dy*dy + dz*dz;
            temp_exp = m_sigma*m_sigma*m_sigma*m_sigma*m_sigma*m_sigma/(r2_ij*r2_ij*r2_ij);
            temp_exp2 = temp_exp*temp_exp;
            m_eps24 = 24*m_epsilon;
            U = 4*m_epsilon*(temp_exp2 - temp_exp);
            m_potentialEnergy += U;
            Fx = m_eps24*(2*temp_exp2 - temp_exp)*dx/(r2_ij);
            Fy = m_eps24*(2*temp_exp2 - temp_exp)*dy/(r2_ij);
            Fz = m_eps24*(2*temp_exp2 - temp_exp)*dz/(r2_ij);
            //std::cout << dx << " " << dy << " " << dz << std::endl;
            atomj->force += {Fx, Fy, Fz}; //This one should be +
            atomi->force -= {Fx, Fy, Fz}; //This one should be - Have checked it with only two atoms
        }
       //atomj->force = {0, 0, 0};
       //std::cout << atomj->force << std::endl;
    }
}
