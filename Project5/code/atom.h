#ifndef ATOM_H
#define ATOM_H
#include "math/vec3.h"

class Atom
{
private:
    float m_mass;   
public:
    vec3 position;
    vec3 velocity;
    vec3 force;
    vec3 initialPosition;
    vec3 m_distanceTravelled;

    Atom(double mass);
    void resetForce();
    void resetVelocityMaxwellian(double temperature);

    double mass() { return m_mass; }
    double m_distanceBeforePBC;
    double m_valueDistanceTravelled;
    double distanceBeforePBC() { return m_distanceBeforePBC; }
    void setMass(double mass) { m_mass = mass; }
};
#endif
