#ifndef PENNINGTRAP_H_
#define PENNINGTRAP_H_

#include "Particle.hpp"
#include <armadillo>
#include <string>
#include <vector>

class PenningTrap{

    public:
        double B0;  // Magnetic field strength
        double V0;  // Voltage
        double d;   // Distance between plates    
        double ke = 1.3893533e5; // Coulomb constant    
        std::vector<Particle> particles;   // Vector of particles
        bool interaction = false;   // Interaction between particles

        // For EX9 not useful for now
        double frequency;   // Frequency of the oscillation
        double amplitude;   // Amplitude of the oscillation

        // Constructor
        PenningTrap(double B0_in, double V0_in, double d_in);


        // Add a particle to the trap 
        //can probably be generalized to add multiple particles at once
        void add_particle(Particle p_in);

        // calculate the external E_field without time dependence
        arma::vec external_E_field(arma::vec r);

        //calculate the external B_field without time dependence
        arma::vec external_B_field(arma::vec r);

        // Force on particle_i from particle_j
        arma::vec force_particle(int i, int j);

        // The total force on particle_i from the external fields Lorentz force
        arma::vec total_force_external(int i);

        // The total force on particle_i from the other particles
        arma::vec total_force_particles(int i);

        // Total force so external force and if interaction is true 
        //particle interaction
        arma::vec total_force(int i);

        // Evolve the system one time step (dt) using Runge-Kutta 4th order
        void evolve_RK4(double dt);

        // Evolve the system one time step (dt) using Forward Euler
        void evolve_forward_Euler(double dt);


};

#endif // PENNINGTRAP_H_