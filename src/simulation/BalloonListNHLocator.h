//  Copyright by Christoph Saffer, Paul Rudolph, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#ifndef BALLOONLISTNHLOCATOR_H
#define    BALLOONLISTNHLOCATOR_H

#include <map>
#include <boost/thread/condition_variable.hpp>

#include "simulation/NeighbourhoodLocator.h"
#include "simulation/Site.h"


class BalloonListNHLocator : public NeighbourhoodLocator {

public:
    /// Class for detecting collisions between agents while maintaining a discrete grid of possible cell positions
    /// Keeps track of the neighbours of a cell throughout a simulation to prevent brute force neighbourhood detection
    BalloonListNHLocator(double gridConstant, Coordinate3D lowerValues, Coordinate3D upperValues, Site *site);

    virtual ~BalloonListNHLocator();

    /// Initializes neighbourhood grid for Site
    void instantiate(Site *site);

    /// Updates neighbourhood lists
    void updateDataStructures(SphereRepresentation *sphereRep) final;

    /// Adds a new sphere to the list
    void addSphereRepresentation(SphereRepresentation *sphereRep) final;

    /// Removes sphere from the list
    void removeSphereRepresentation(SphereRepresentation *sphereRep) final;

    int controlFunction() final;
    bool hasCollision(Agent *agent) final;
    std::vector<std::shared_ptr<Collision>> getCollisions(Agent *agent) final;
    std::vector<Coordinate3D> getCollisionSpheres(SphereRepresentation *sphereRep, Coordinate3D dirVec) final;
    std::string getTypeName() final;
    int getNumberOfAgentTypeInBalloonList(std::string agentType) final;

private:
    double gridConstant;
    std::vector<std::vector<std::vector<std::vector<SphereRepresentation *> > > > balloonList;
    std::map<SphereRepresentation *, std::vector<int> > sphereRepresentationAllocator;
    Coordinate3D lowerPoint;
    Coordinate3D upperPoint;
    int gridSize[3];
    int complexes;
    double timestep;
    boost::condition_variable m_cond;
    int checksum;
    void initialGridCreation();
    bool checkCollisions(std::vector<std::shared_ptr<Collision> > *neighbours,Agent *agent,SphereRepresentation *sphereRep,
                         int u, int v, int w, bool justCheck = false);
};

#endif    /* BALLOONLISTNHLOCATOR_H */