#ifndef PARTICLE_H_
#define PARTICLE_H_

#include <armadillo>
#include <string>
#include <vector>

// Structure
struct Particle{
    Particle(double q, double m, arma::vec r, arma::vec v);
    double charge;
    double mass;
    arma::vec position;
    arma::vec velocity;

    

};

#endif // PARTICLE_H_