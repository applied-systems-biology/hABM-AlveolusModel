//  Copyright by Christoph Saffer, Paul Rudolph, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#include "simulation/Particle.h"
#include "simulation/Site.h"
#include "simulation/movement/RandomWalk.h"
#include "utils/macros.h"


Particle::Particle(unsigned int id, Site *site, const Coordinate3D &pos, double conc, double dc) {

    this->id = id;
    this->position = std::make_unique<Coordinate3D>();
    this->initPosition = pos;
    *this->position = pos;
    this->concentration = conc;
    this->dc = dc;
    this->concChangeInCurTimestep = 0;
    this->movement = std::make_unique<RandomWalk>(site, site->getNumberOfSpatialDimensions());
    movement->setSite(site);
    movement->setCurrentPosition(position.get());
    particleNeighbourList = std::make_unique<ParticleNeighbourList>(this);
    isInSite = site->containsPosition(*position);
}

void Particle::doAllActionsForTimestep(double timestep) {
    diffusePSE(timestep);
}

double Particle::estimateLowestTimestep() {
    std::vector<Particle *> nParticle = getParticleNeighbourList()->getNeighbours();
    std::vector<double> distances = getParticleNeighbourList()->getDistances();
    std::vector<double> contactAreas = getParticleNeighbourList()->getContactAreas();

    unsigned int size = distances.size();

    double minTimestep = 1000000;
    double sum = 0;
    for (unsigned int i = 0; i < size; i++) {
        sum += (dc * contactAreas[i]) / (distances[i] * area);
    }
    minTimestep = 1.0 / sum;
    return minTimestep;
}

void Particle::diffusePSE(double timestep) {
    if (isInSite) { //only "in site" grid points are used for the calculations
        std::vector<Particle *> neighbours = particleNeighbourList->getNeighbours();
        std::vector<double> preFactorsPSE = particleNeighbourList->getPreFactorsPSE();
        std::vector<Particle *>::iterator it = neighbours.begin();
        std::vector<double>::iterator itPF = preFactorsPSE.begin();
        double curOwnPrefactor = 0;
        double concChangeDiffusion = 0;
        while (it != neighbours.end()) {
            concChangeDiffusion += (*it)->getConcentration() * (*itPF);
            curOwnPrefactor += (*itPF);
            it++;
            itPF++;
        }
        concChangeDiffusion -= concentration * curOwnPrefactor;
        concChangeDiffusion *= timestep;

        addConcentrationChange(concChangeDiffusion);
    }
}

Coordinate3D Particle::getGradient() {
    Coordinate3D gradient = Coordinate3D();
    if (isInSite) { //only "in site" grid points are used for the calculations
        std::vector<Particle *> neighbours = particleNeighbourList->getNeighbours();
        std::vector<double> preFactorsGradient = particleNeighbourList->getPreFactorsGradient();
        std::vector<Particle *>::iterator it = neighbours.begin();
        unsigned int index = 0;
        while (it != neighbours.end()) {

            //Sukumar et al.
            Coordinate3D curGradient{(*it)->getPosition() - *position};
            curGradient *= 0.5 * preFactorsGradient[index] * ((*it)->getConcentration() - concentration);
            it++;
            index++;
            gradient += curGradient;
        }
    }

    return gradient;
}

double Particle::getGradientStrength() {
    Coordinate3D gradient = Coordinate3D();
    if (isInSite) {
        std::vector<Particle *> neighbours = particleNeighbourList->getNeighbours();
        std::vector<double> preFactorsGradient = particleNeighbourList->getPreFactorsGradient();

        std::vector<Particle *>::iterator it = neighbours.begin();
        unsigned int index = 0;
        while (it != neighbours.end()) {

            //Sukumar et al.
            Coordinate3D curGradient{(*it)->getPosition() - *position};
            curGradient *= 0.5 * preFactorsGradient[index] * ((*it)->getConcentration() - concentration);
            it++;
            index++;
            gradient += curGradient;
        }
    }

    return gradient.getMagnitude();
}

double Particle::getGradientDirection() {
    Coordinate3D gradient = Coordinate3D();

    if (isInSite) {
        std::vector<Particle *> neighbours = particleNeighbourList->getNeighbours();
        std::vector<double> preFactorsGradient = particleNeighbourList->getPreFactorsGradient();
        std::vector<Particle *>::iterator it = neighbours.begin();
        unsigned int index = 0;
        while (it != neighbours.end()) {

            //Sukumar et al.
            Coordinate3D curGradient{(*it)->getPosition() - *position};
            curGradient *= 0.5 * preFactorsGradient[index] * ((*it)->getConcentration() - concentration);
            it++;
            index++;
            gradient += curGradient;
        }
    }

    gradient.setMagnitude(1.0);
    Coordinate3D gradientPos = getPosition();
    gradientPos += gradient;

    SphericCoordinate3D ownPos = abm::util::toSphericCoordinates(*position);
    SphericCoordinate3D goalPos = abm::util::toSphericCoordinates(gradientPos);

    //get the alpha turning angle by taking the spheric coordinates as input:
    double alpha = 0;
    double dphi = ownPos.phi - goalPos.phi;
    double dphiAbs = fabs(dphi);
    double c = ownPos.theta;
    double b = goalPos.theta;
    double a = acos(sin(ownPos.theta) * sin(goalPos.theta) * cos(dphi) + cos(ownPos.theta) * cos(goalPos.theta));

    if (c == 0) {
        alpha = goalPos.phi;
    } else {
        if (b == 0) {
            alpha = M_PI;
        } else {
            //compute beta by cases
            double beta;
            if (ownPos.phi == goalPos.phi) {
                if (ownPos.theta > goalPos.theta) {
                    beta = 0;
                } else {
                    beta = M_PI;
                }
            } else {
                beta = acos((cos(b) - cos(a) * cos(c)) / (sin(a) * sin(c)));
            }

            //compute alpha by cases
            if ((ownPos.phi <= goalPos.phi && dphiAbs <= M_PI) || (ownPos.phi > goalPos.phi && dphiAbs > M_PI)) {
                alpha = M_PI - beta;
            } else {
                alpha = M_PI + beta;
            }

        }
    }

    return alpha;
}

double Particle::applyConcentrationChange(double time_delta) {

    double relChange = 0;
    if (isInSite) { //only "in site" grid points are used for the calculations
        if (concentration > 0) {
            relChange = fabs(concChangeInCurTimestep / (concentration * time_delta));
        }
        concentration += concChangeInCurTimestep;
        concChangeInCurTimestep = 0;
    }
    return relChange;
}

void Particle::setLabel() {
    if (isInSite) {
        isBoundaryParticle = false;
    } else {
        std::vector<Particle *> neighbours = particleNeighbourList->getNeighbours();
        std::vector<Particle *>::iterator itN = neighbours.begin();

        isBoundaryParticle = false;

        while (itN != neighbours.end()) {
            if ((*itN)->isInSite) {
                isBoundaryParticle = true;
                break;
            }
            itN++;
        }
    }
}