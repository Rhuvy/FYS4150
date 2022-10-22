#include "Particle.hpp"


Particle::Particle(double q, double m, arma::vec r, arma::vec v){
    charge = q;
    mass = m;
    position = r;
    velocity = v;
}

double Particle::get_charge(){
    return charge;
}

double Particle::get_mass(){
    return mass;
}


arma::vec Particle::get_position(){
    return position;
}


arma::vec Particle::get_velocity(){
    return velocity;
}


void Particle::set_position(arma::vec r){
    position = r;
}


void Particle::set_velocity(arma::vec v){
    velocity = v;
}

