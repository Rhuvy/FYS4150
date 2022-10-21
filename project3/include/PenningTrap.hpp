#ifndef PENNINGTRAP_H_
#define PENNINGTRAP_H_

#include "Particle.hpp"
#include <armadillo>
#include <string>
#include <vector>

struct PenningTrap{
    PenningTrap(std::vector<Particle> particles, double B0_in, double V0_in, double d_in);
    double B0_;
    double V0_;
    double d_;
    int N;
    bool interaction;
    std::vector<Particle> particles_;



};

#endif // PENNINGTRAP_H_