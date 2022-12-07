//  Copyright by Christoph Saffer, Paul Rudolph, Sandra Timme, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#ifndef NEIGHBOURHOODLOCATOR_H
#define    NEIGHBOURHOODLOCATOR_H

#include <set>
#include <algorithm>
#include <string>
#include <vector>

#include "simulation/morphology/SphereRepresentation.h"


class Agent;
class Collision;
class Site;

class NeighbourhoodLocator {
public:
    /// Abstract class for a cell neighbourhood locator to speed up collision checks between cells
    NeighbourhoodLocator(Site *Site);
    virtual ~NeighbourhoodLocator();
    virtual void instantiate();
    virtual std::vector<std::shared_ptr<Collision>> getCollisions(Agent *agent);
    virtual bool hasCollision(Agent *agent) { return false; };
    virtual std::vector<Coordinate3D> getCollisionSpheres(SphereRepresentation *sphereRep, Coordinate3D dirVec);
    virtual void updateDataStructures(SphereRepresentation *sphereRep);
    virtual void setThresholdDistance(double thresh);
    void setInteractionCheckInterval(unsigned int icInterval) { checkInteractionsTimestepInterval = icInterval; }
    unsigned int getInteractionCheckInterval() { return checkInteractionsTimestepInterval; };
    virtual void removeSphereRepresentation(SphereRepresentation *sphereRep);
    virtual void addSphereRepresentation(SphereRepresentation *sphereRep);
    virtual int controlFunction() { return 0; };
    virtual std::string getTypeName();
    virtual void check() {};
    virtual int getNumberOfAgentTypeInBalloonList(std::string agentType) { return 0; };
protected:
    double thresholdDistance;
    unsigned int checkInteractionsTimestepInterval;
    Site *site_;

};

#endif    /* NEIGHBOURHOODLOCATOR_H */

