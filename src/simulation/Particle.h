//  Copyright by Christoph Saffer, Paul Rudolph, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#ifndef PARTICLE_H
#define    PARTICLE_H

#include <memory>

#include "basic/Coordinate3D.h"
#include "basic/SphericCoordinate3D.h"
#include "simulation/ParticleNeighbourList.h"

class Site;
class RandomWalk;

class Particle {

public:
    /// Class for modelling particles in a 'Lagrangian' way
    Particle() = default;
    Particle(unsigned int id, Site *site, const Coordinate3D &pos, double conc, double dc);

    /*!
     * Performs all actions for one timestep for one particle
     * @param timestep Double that contains timestep
     */
    void doAllActionsForTimestep(double timestep);

    /*!
     * Performs Particle Strength Exchange (PSE) or concentration diffusion on single particle
     * @param timestep
     */
    void diffusePSE(double timestep);

    /*!
     * Adds concentration change to particle
     * @param value Double that contains value that is added to concentration of particle
     */
    void addConcentrationChange(double value) { concChangeInCurTimestep += value; }

    /*!
     * Adds and maintains neighbour list for one particle
     * @param p Particle that contains particle that is added to neighbour list
     * @param contactArea Double that contains contacting area
     */
    void addNeighbour(Particle *p, double contactArea) { particleNeighbourList->addNeighbour(p, contactArea); };

    /*!
     * Applys concentration change for timestep
     * @param time_delta Double that contains timestep
     * @return Double that concentration change
     */
    double applyConcentrationChange(double time_delta);

    /*!
     * Calculates lowest timestep possible for current configuration
     * @return Double that contains timestep
     */
    double estimateLowestTimestep();

    void setArea(double areaP) { area = areaP; };
    void setLabel();
    double getGradientStrength();
    double getGradientDirection();
    double getDiffusionCoefficient() const { return dc; };
    double *getConcentrationRef() { return &concentration; };
    Coordinate3D getPosition() { return *position; };
    Coordinate3D getGradient();
    ParticleNeighbourList *getParticleNeighbourList() { return particleNeighbourList.get(); }
    [[nodiscard]] unsigned int getId() const { return id; };
    [[nodiscard]] double getArea() const { return area; };
    [[nodiscard]] double getConcentration() const { return concentration; };
    [[nodiscard]] bool getIsInSite() const { return isInSite; };
    [[nodiscard]] bool getIsAtBoundary() const { return isBoundaryParticle; };

private:
    unsigned int id;
    Site *site{};
    bool isStatic{};
    bool isInSite;
    bool isBoundaryParticle{};
    double area{};
    std::unique_ptr<Coordinate3D> position;
    Coordinate3D initPosition;
    double concentration;
    double concChangeInCurTimestep;
    std::unique_ptr<RandomWalk> movement;
    double dc;
    std::unique_ptr<ParticleNeighbourList> particleNeighbourList;

};

#endif    /* PARTICLE_H */

