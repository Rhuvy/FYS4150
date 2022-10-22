#include "PenningTrap.hpp"
#include "Particle.hpp" 
  
// Constructor
PenningTrap::PenningTrap(double B0_in, double V0_in, double d_in){
  B0 = B0_in;
  V0 = V0_in;
  d = d_in;
}

// Add a particle to the trap
void PenningTrap::add_particle(Particle p_in){
  PenningTrap::particles.push_back(p_in);
}

  // External electric field at point r=(x,y,z)
arma::vec PenningTrap::external_E_field(arma::vec r){
  arma::vec E = {0, 0, 0};
    
  if(arma::norm(r) < d){

    // Ugly but more readable if there is a problem 
    // could define V/d^2 but just want to make a code working 
    E(0) = r(0) * V0 / (d * d);
    E(1) = r(1) * V0 / (d * d);
    E(2) = -2. * r(2) * V0 / (d * d);
  }

  return E;
}

// External magnetic field at point r=(x,y,z)
arma::vec PenningTrap::external_B_field(arma::vec r){
  
  // WE check the norm because we dont have B outside the trap
  if (arma::norm(r) < d){
    return {0, 0, B0};
  }
 
  else{
    return {0, 0, 0};
  }
}

// Force on particle_i from particle_j
arma::vec PenningTrap::force_particle(int i, int j){
  
  // Calculate the distance between the particles and its norm
  arma::vec r_ij = particles[i].get_position() - particles[j].get_position();
  double r = arma::norm(r_ij);

  // The charge of each particle
  double q_i = particles[i].get_charge();
  double q_j = particles[j].get_charge();
  
  // Initialize the force 
  arma::vec F = {0, 0, 0};
  
  // Check if the particles are not the same maybe it is the problem
  if (i != j){
    //Ugly formula but more readable because i dont know where the fuck our problem is 
    F = ke * (q_i * q_j * r_ij) / (r * r * r);
  }
  
  return F;

}


// The total force on particle_i from the external fields
arma::vec PenningTrap::total_force_external(int i){
  
  arma::vec F = {0, 0, 0};
  
  // Calculate the force on particle i from the external fields
  F = particles[i].get_charge() * (external_E_field(particles[i].get_position()) + arma::cross(particles[i].get_velocity(), external_B_field(particles[i].get_position())));
  
  return F;

}

// The total force on particle_i from the other particles
arma::vec PenningTrap::total_force_particles(int i){
    
    arma::vec F = {0, 0, 0};
    
    // Calculate the force on particle i from the other particles
    for (int j = 0; j < particles.size(); j++){
      F += force_particle(i, j);
    }
    
    return F;
  
}

// The total force on particle_i from both external fields and other particles
arma::vec PenningTrap::total_force(int i){

  arma::vec F = {0, 0, 0};
  // Calculate the total force on particle i
  F = total_force_external(i);

  if(interaction);
    F += total_force_particles(i);

  return F;
}

// Evolve the system one time step (dt) using Runge-Kutta 4th order
void PenningTrap::evolve_RK4(double dt){
    // IMPORTANT
    // We have to loop separetly the coefficients 
    // because we need to use the old positions and velocities to calculate the new ones
    // when we have interaction since we need the old position at all time 
    // in the calculation of the force


    // Initialization
    int N = particles.size();

    // Create a initial 'pictures' of our particle 
    // Like this we can always access the intial position and velocity
    std::vector<Particle> particles_initial = particles;
 
    // k is for the velocity and l for the position
    std::vector<arma::vec> k1(N), k2(N), k3(N), k4(N);
    std::vector<arma::vec> l1(N), l2(N), l3(N), l4(N);

    // Calculate the k1 and l1 for each particle
    for (int i = 0; i < N; i++){
      // Calculate k1 and l1
      k1[i] = dt * total_force(i) / particles[i].get_mass();
      l1[i] = dt * particles[i].get_velocity();
    }

    for (int i = 0; i < N; i++){
      //Calculate new r and v
      arma::vec r = particles_initial[i].get_position() + l1[i] / 2;
      arma::vec v = particles_initial[i].get_velocity() + k1[i] / 2;

      // Update the particle
      particles[i].set_position(r);
      particles[i].set_velocity(v);
    }

    // Calculate the k2 and l2 for each particle
    for (int i = 0; i < N; i++){
      // Calculate k2 and l2
      k2[i] = dt * total_force(i) / particles[i].get_mass();
      l2[i] = dt * particles[i].get_velocity();
    }

    for (int i = 0; i < N; i++){
      //Calculate new r and v
      arma::vec r = particles_initial[i].get_position() + l2[i] / 2;
      arma::vec v = particles_initial[i].get_velocity() + k2[i] / 2;

      // Update the particle
      particles[i].set_position(r);
      particles[i].set_velocity(v);
    }

    // Calculate the k3 and l3 for each particle
    for (int i = 0; i < N; i++){
      // Calculate k3 and l3
      k3[i] = dt * total_force(i) / particles[i].get_mass();
      l3[i] = dt * particles[i].get_velocity();
    }

    for (int i = 0; i < N; i++){
      //Calculate new r and v
      arma::vec r = particles_initial[i].get_position() + l3[i];
      arma::vec v = particles_initial[i].get_velocity() + k3[i];

      // Update the particle
      particles[i].set_position(r);
      particles[i].set_velocity(v);
    }

    // Calculate the k4 and l4 for each particle

    for (int i = 0; i < N; i++){
      // Calculate k4 and l4
      k4[i] = dt * total_force(i) / particles[i].get_mass();
      l4[i] = dt * particles[i].get_velocity();
    }

    for (int i = 0; i < N; i++){

      //Calculate new r and v
      arma::vec r = particles_initial[i].get_position() + (l1[i] + 2 * l2[i] + 2 * l3[i] + l4[i]) / 6;
      arma::vec v = particles_initial[i].get_velocity() + (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6;

      // Update the particle
      particles[i].set_position(r);
      particles[i].set_velocity(v);
    }

}

// Evolve the system one time step (dt) using Forward Euler
void PenningTrap::evolve_forward_Euler(double dt){
  
  // Initialization
  int N = particles.size();
  arma::mat r_new = arma::zeros<arma::mat>(N, 3);
  arma::mat v_new = arma::zeros<arma::mat>(N, 3);

  for(int i = 0; i < N; i++){

    //Current position and speed of particle i
    arma::vec r_i = particles[i].get_position();
    arma::vec v_i = particles[i].get_velocity();

    // Calculate the acceleration so F/m
    arma::vec F = total_force(i);
    arma::vec a = F / particles[i].get_mass();

    // Calculate the new position and velocity of each particle
    r_new.row(i) = r_i + dt * v_i;
    v_new.row(i) = v_i + dt * a;
  }

  //Now we update the position and velocity of each particle in particles
  for(int i = 0; i < N; i++){
    particles[i].set_position(r_new.row(i));
    particles[i].set_velocity(v_new.row(i));
  }

}

