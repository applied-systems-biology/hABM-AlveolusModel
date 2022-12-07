//  Copyright by Christoph Saffer, Paul Rudolph, Sandra Timme, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#ifndef PARTICLENEIGHBOURLIST_H
#define    PARTICLENEIGHBOURLIST_H

#include <vector>

class Particle;

class ParticleNeighbourList {
  public:
    /// Class for implementation of grid-based particle interaction and diffusion
    ParticleNeighbourList() {};
    ParticleNeighbourList(Particle *p) { particle = p; };

    void addNeighbour(Particle *particle, double contactArea);
    bool existsInList(Particle *p);

    std::vector<Particle *> &getNeighbours() { return neighbourParticles; };
    std::vector<double> &getDistances() { return neighbourDistances; };
    std::vector<double> &getContactAreas() { return neighbourContactArea; };
    std::vector<double> &getPreFactorsPSE() { return preFactorsPSE; };
    std::vector<double> &getPreFactorsGradient() { return preFactorsGradient; };

private:
    Particle *particle;
    std::vector<Particle *> neighbourParticles;
    std::vector<double> neighbourDistances;
    std::vector<double> neighbourContactArea;
    std::vector<double> preFactorsPSE;
    std::vector<double> preFactorsGradient;
};

#endif    /* PARTICLENEIGHBOURLIST_H */

