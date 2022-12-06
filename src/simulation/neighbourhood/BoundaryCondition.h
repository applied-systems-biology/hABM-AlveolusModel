//  Copyright by Christoph Saffer, Paul Rudolph, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#ifndef BOUNDARYCONDITION_H
#define    BOUNDARYCONDITION_H

#include <memory>

#include "simulation/Agent.h"
#include "basic/Coordinate3D.h"
#include "simulation/Particle.h"

class Analyser;
class Site;

class BoundaryCondition {
public:
  // Abstract class for boundary condition
    explicit BoundaryCondition(Site *site);

    virtual ~BoundaryCondition() = default;
    virtual void handleBoundaryCross(Agent *agent, Coordinate3D *moveVec, double current_time);
    virtual void handleBoundaryCross(Particle *particle, Coordinate3D *moveVec) {};
    virtual std::string getTypeName();

protected:
    Site *site;
};

#endif    /* BOUNDARYCONDITION_H */

