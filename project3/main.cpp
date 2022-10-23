#include <iostream>
#include <armadillo>
#include "PenningTrap.hpp"
#include "Particle.hpp"

////// A mess just to check everything working needs some polishing but normally everythings works as planned //////

int main(){

    // Initialize the particle
    double q = 1;
    double m = 40.08;

    arma::vec r_initial = {20., 0., 20.};
    arma::vec v_initial = {0., 25., 0.};

    arma::vec r2_initial = {25., 25., 0.};
    arma::vec v2_initial = {0., 40., 5.};

    Particle p1(q, m, r_initial, v_initial);
    Particle p2(q, m, r2_initial, v2_initial);

    // Initialize the trap
    double B0 = 9.65e1;
    double V0 = 2.41e6;
    double d = 500.;

    PenningTrap trap(B0, V0, d);
    // trap.interaction = true;
    // add the particle to the trap 
    trap.add_particle(p1);
    trap.add_particle(p2);

    // Initialize the time and position 
    double dt = 0.001;
    int N = 50000; // should be 50 ms / dt to test

    arma::vec t = arma::linspace(0, N*dt, N);
    arma::vec x = arma::vec(N, arma::fill::zeros);
    arma::vec y = arma::vec(N, arma::fill::zeros);
    arma::vec z = arma::vec(N, arma::fill::zeros);

    arma::vec v_x = arma::vec(N, arma::fill::zeros);
    arma::vec v_y = arma::vec(N, arma::fill::zeros);
    arma::vec v_z = arma::vec(N, arma::fill::zeros);

    arma::vec x2 = arma::vec(N, arma::fill::zeros);
    arma::vec y2 = arma::vec(N, arma::fill::zeros);
    arma::vec z2 = arma::vec(N, arma::fill::zeros);

    arma::vec v_x2 = arma::vec(N, arma::fill::zeros);
    arma::vec v_y2 = arma::vec(N, arma::fill::zeros);
    arma::vec v_z2 = arma::vec(N, arma::fill::zeros);

    // Initialize the value 
    x(0) = p1.get_position()(0);
    y(0) = p1.get_position()(1);
    z(0) = p1.get_position()(2);

    v_x(0) = p1.get_velocity()(0);
    v_y(0) = p1.get_velocity()(1);
    v_z(0) = p1.get_velocity()(2);

    x2(0) = p2.get_position()(0);
    y2(0) = p2.get_position()(1);
    z2(0) = p2.get_position()(2);

    v_x2(0) = p2.get_velocity()(0);
    v_y2(0) = p2.get_velocity()(1);
    v_z2(0) = p2.get_velocity()(2);

    // For RK4
    // arma::vec x4 = arma::vec(N, arma::fill::zeros);
    // arma::vec y4 = arma::vec(N, arma::fill::zeros);
    // arma::vec z4 = arma::vec(N, arma::fill::zeros);

    // Initialize the value
    // x4(0) = p1.get_position()(0);
    // y4(0) = p1.get_position()(1);
    // z4(0) = p1.get_position()(2);

    // Test FE 
    for (int i = 0; i < N-1; i++){
        
        // move one step forward
        trap.evolve_RK4(dt);

        // update the position
        x(i+1) = trap.particles[0].get_position()(0);
        y(i+1) = trap.particles[0].get_position()(1);
        z(i+1) = trap.particles[0].get_position()(2);

        v_x(i+1) = trap.particles[0].get_velocity()(0);
        v_y(i+1) = trap.particles[0].get_velocity()(1);
        v_z(i+1) = trap.particles[0].get_velocity()(2);

        x2(i+1) = trap.particles[1].get_position()(0);
        y2(i+1) = trap.particles[1].get_position()(1);
        z2(i+1) = trap.particles[1].get_position()(2);

        v_x2(i+1) = trap.particles[1].get_velocity()(0);
        v_y2(i+1) = trap.particles[1].get_velocity()(1);
        v_z2(i+1) = trap.particles[1].get_velocity()(2);
    }

    // Write to file
    x.save("x_RK4_no.txt", arma::raw_ascii);
    y.save("y_RK4_no.txt", arma::raw_ascii);
    z.save("z_RK4_no.txt", arma::raw_ascii);

    v_x.save("v_x_RK4_no.txt", arma::raw_ascii);
    v_y.save("v_y_RK4_no.txt", arma::raw_ascii);
    v_z.save("v_z_RK4_no.txt", arma::raw_ascii);

    x2.save("x2_RK4_no.txt", arma::raw_ascii);
    y2.save("y2_RK4_no.txt", arma::raw_ascii);
    z2.save("z2_RK4_no.txt", arma::raw_ascii);

    v_x2.save("v_x2_RK4_no.txt", arma::raw_ascii);
    v_y2.save("v_y2_RK4_no.txt", arma::raw_ascii);
    v_z2.save("v_z2_RK4_no.txt", arma::raw_ascii);


    // // Test RK4
    // PenningTrap trap2(B0, V0, d);
    // trap2.add_particle(p1);

    // for (int i = 0; i < N-1; i++){

    //     // move one step forward
    //     trap2.evolve_RK4(dt);

    //     // update the position
    //     x4(i+1) = trap2.particles[0].get_position()(0);
    //     y4(i+1) = trap2.particles[0].get_position()(1);
    //     z4(i+1) = trap2.particles[0].get_position()(2);
    // }

    // // Write to file
    // x4.save("x_RK4.txt", arma::raw_ascii);
    // y4.save("y_RK4.txt", arma::raw_ascii);
    // z4.save("z_RK4.txt", arma::raw_ascii);

    // Two particle system with and without interaction
    // PenningTrap trap3(B0, V0, d);
    // trap3.add_particle(p1);
    // trap3.add_particle(p2);

    // arma::vec x3_1 = arma::vec(N, arma::fill::zeros);
    // arma::vec y3_1 = arma::vec(N, arma::fill::zeros);
    // arma::vec z3_1 = arma::vec(N, arma::fill::zeros);

    // arma::vec vx3_1 = arma::vec(N, arma::fill::zeros);
    // arma::vec vy3_1 = arma::vec(N, arma::fill::zeros);
    // arma::vec vz3_1 = arma::vec(N, arma::fill::zeros);

    // arma::vec vx3_2 = arma::vec(N, arma::fill::zeros);
    // arma::vec vy3_2 = arma::vec(N, arma::fill::zeros);
    // arma::vec vz3_2 = arma::vec(N, arma::fill::zeros);

    // arma::vec x3_2 = arma::vec(N, arma::fill::zeros);
    // arma::vec y3_2 = arma::vec(N, arma::fill::zeros);
    // arma::vec z3_2 = arma::vec(N, arma::fill::zeros);


    return 0;
}
