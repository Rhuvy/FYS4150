#ifndef PARTICLE_H_
#define PARTICLE_H_

#include <armadillo>
#include <string>
#include <vector>

// Structure
class Particle{

    // Definition of public so that we can access the variables from outside the class (?)
    public:
        double charge;
        double mass;
        arma::vec position;
        arma::vec velocity;

    Particle(double q, double m, arma::vec r, arma::vec v);

    double get_charge();    // Returns the charge of the particle
    
    double get_mass();    // Returns the mass of the particle 

    arma::vec get_position();   // Returns the position of the particle

    arma::vec get_velocity();   // Returns the velocity of the particle

    void set_position(arma::vec r);   // Sets the position of the particle

    void set_velocity(arma::vec v);   // Sets the velocity of the particle


};

#endif // PARTICLE_H_