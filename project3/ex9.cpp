#include "PenningTrap.hpp"
#include "Particle.hpp"
#include <iostream>
#include<armadillo>

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
    trap.frequency = 0.2;
    trap.amplitude = 0;

    // this function works 

    // trap.add_particles(n);
    // double test = trap.number_of_particles();
    // std::cout << test << std::endl;

    // trap.add_particles(n);
    double time = 50.;
    double dt = 0.001;
    // Initialize the time and position
    int N = time/dt;

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

    // Runge-Kutta 4th order
    for (int i = 0; i < N-1; i++){
        trap.evolve_RK4(dt, t(i));

        // Update the position and velocity
        x(i+1) = p1.get_position()(0);
        y(i+1) = p1.get_position()(1);
        z(i+1) = p1.get_position()(2);

        v_x(i+1) = p1.get_velocity()(0);
        v_y(i+1) = p1.get_velocity()(1);
        v_z(i+1) = p1.get_velocity()(2);

        x2(i+1) = p2.get_position()(0);
        y2(i+1) = p2.get_position()(1);
        z2(i+1) = p2.get_position()(2);

        v_x2(i+1) = p2.get_velocity()(0);
        v_y2(i+1) = p2.get_velocity()(1);
        v_z2(i+1) = p2.get_velocity()(2);

    }
    //Write to file
    x.save("xf.txt", arma::raw_ascii);
    y.save("yf.txt", arma::raw_ascii);
    z.save("zf.txt", arma::raw_ascii);

    v_x.save("vx_f.txt", arma::raw_ascii);
    v_y.save("vy_f.txt", arma::raw_ascii);
    v_z.save("vz_f.txt", arma::raw_ascii);

    x2.save("x2_f.txt", arma::raw_ascii);
    y2.save("y2_f.txt", arma::raw_ascii);
    z2.save("z2_f.txt", arma::raw_ascii);

    v_x2.save("vx2_f.txt", arma::raw_ascii);
    v_y2.save("vy2_f.txt", arma::raw_ascii);
    v_z2.save("vz2_f.txt", arma::raw_ascii);

    return 0;
}
