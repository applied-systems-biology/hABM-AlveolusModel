//  Copyright by Christoph Saffer, Paul Rudolph, Sandra Timme, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#include <vector>

#include "simulation/ParticleNeighbourList.h"
#include "simulation/Particle.h"


void ParticleNeighbourList::addNeighbour(Particle *p, double contactArea) {

    double distance;
    double preFactorPSE;
    double preFactorGradient;

    if (!existsInList(p)) {
        neighbourParticles.push_back(p);
        distance = p->getPosition().calculateEuclidianDistance(particle->getPosition());
        neighbourDistances.push_back(distance);
        neighbourContactArea.push_back(contactArea);

        preFactorPSE = particle->getDiffusionCoefficient() * contactArea / (distance * particle->getArea());
        preFactorsPSE.push_back(preFactorPSE);

        preFactorGradient = contactArea / (distance * particle->getArea());
        preFactorsGradient.push_back(preFactorGradient);

    }
}

bool ParticleNeighbourList::existsInList(Particle *p) {
    bool isIn = false;

    auto it = neighbourParticles.begin();
    while (it != neighbourParticles.end()) {
        if ((*it) == p) {
            isIn = true;
            break;
        }
        it++;
    }
    return isIn;
}