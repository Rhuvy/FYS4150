#include <iostream>
#include <armadillo>
#include "PenningTrap.hpp"
#include "Particle.hpp"


int main(){

    // Initialize the particle
    double q = 1;
    double m = 40.08;

    arma::vec r_initial = {20., 0., 20.};
    arma::vec v_initial = {0., 25., 0.};

    Particle p1(q, m, r_initial, v_initial);

    // Initialize the trap
    double B0 = 9.65e1;
    double V0 = 2.41e6;
    double d = 500.;

    PenningTrap trap(B0, V0, d);
    
    // add the particle to the trap 
    trap.add_particle(p1);

    // Initialize the time and position 
    double dt = 0.001;
    int N = 50000; // should be 50 ms / dt to test

    arma::vec t = arma::linspace(0, N*dt, N);
    arma::vec x = arma::zeros<arma::vec>(N);
    arma::vec y = arma::zeros<arma::vec>(N);
    arma::vec z = arma::zeros<arma::vec>(N);

    // Initialize the value 
    x(0) = p1.get_position()(0);
    y(0) = p1.get_position()(1);
    z(0) = p1.get_position()(2);

    // For RK4
    arma::vec x4 = arma::zeros<arma::vec>(N);
    arma::vec y4 = arma::zeros<arma::vec>(N);
    arma::vec z4 = arma::zeros<arma::vec>(N);

    // Initialize the value
    x4(0) = p1.get_position()(0);
    y4(0) = p1.get_position()(1);
    z4(0) = p1.get_position()(2);

    // Test FE 
    for (int i = 0; i < N-1; i++){

        // move one step forward
        trap.evolve_forward_Euler(dt);

        // update the position
        x(i+1) = trap.particles[0].get_position()(0);
        y(i+1) = trap.particles[0].get_position()(1);
        z(i+1) = trap.particles[0].get_position()(2);
    }

    // Write to file
    x.save("x_FE.txt", arma::raw_ascii);
    y.save("y_FE.txt", arma::raw_ascii);
    z.save("z_FE.txt", arma::raw_ascii);

    // Test RK4
    PenningTrap trap2(B0, V0, d);
    trap2.add_particle(p1);

    for (int i = 0; i < N-1; i++){

        // move one step forward
        trap2.evolve_RK4(dt);

        // update the position
        x4(i+1) = trap2.particles[0].get_position()(0);
        y4(i+1) = trap2.particles[0].get_position()(1);
        z4(i+1) = trap2.particles[0].get_position()(2);
    }

    // Write to file
    x4.save("x_RK4.txt", arma::raw_ascii);
    y4.save("y_RK4.txt", arma::raw_ascii);
    z4.save("z_RK4.txt", arma::raw_ascii);


    return 0;
}
