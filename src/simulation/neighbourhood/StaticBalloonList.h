//  Copyright by Christoph Saffer, Paul Rudolph, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#ifndef STATICBALLOONLIST_H
#define    STATICBALLOONLIST_H

#include <vector>
#include <map>

#include "basic/Coordinate3D.h"

class StaticBalloonList {
public:
  // Class for defining a baloon list as a grid that represents the whole environment for efficient neighbourhood detection of cells.
    StaticBalloonList();
    StaticBalloonList(const StaticBalloonList &orig);
    virtual ~StaticBalloonList();

    StaticBalloonList(double gridConstant, Coordinate3D lowerValues, Coordinate3D upperValues);
    void instantiate();
    void getInteractions(Coordinate3D myPos, std::vector<unsigned int> &neighbours);
    unsigned int getClosestObjectIndex(Coordinate3D myPos);
    void getClosestObjectIndices(Coordinate3D myPos, std::vector<unsigned int> &neighbourList,
                                 unsigned int closestXParticles);
    void addCoordinateWithId(Coordinate3D input, unsigned int id);
    void setThreshold(double thresh) { threshold = thresh; };

private:
    double gridConstant;
    double threshold;
    std::vector<std::vector<std::vector<std::vector<unsigned int> > > > balloonList;
    std::map<unsigned int, int *> objectAllocator;
    std::map<unsigned int, Coordinate3D> coordinateAllocator;
    Coordinate3D lowerPoint;
    Coordinate3D upperPoint;
    int gridSize[3];
    void initialGridCreation();
};

#endif    /* STATICBALLOONLIST_H */

