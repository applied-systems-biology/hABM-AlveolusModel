//  Copyright by Christoph Saffer, Paul Rudolph, Sandra Timme, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#include <algorithm>
#include <stdlib.h>
#include <cmath>

#include "StaticBalloonList.h"
#include "utils/macros.h"

StaticBalloonList::StaticBalloonList() {
}

StaticBalloonList::StaticBalloonList(const StaticBalloonList &orig) {
}

StaticBalloonList::~StaticBalloonList() {
}

StaticBalloonList::StaticBalloonList(double gridConstant, Coordinate3D lowerValues, Coordinate3D upperValues) {
    this->gridConstant = gridConstant;
    this->lowerPoint = lowerValues;
    this->upperPoint = upperValues;

    instantiate();
}

void StaticBalloonList::instantiate() {
    initialGridCreation();
}

void StaticBalloonList::initialGridCreation() {
    unsigned int nx, ny, nz;
    nx = (unsigned int) ceil((upperPoint.x - lowerPoint.x) / gridConstant) + 1;
    ny = (unsigned int) ceil((upperPoint.y - lowerPoint.y) / gridConstant) + 1;
    nz = (unsigned int) ceil((upperPoint.z - lowerPoint.z) / gridConstant) + 1;

    gridSize[0] = nx;
    gridSize[1] = ny;
    gridSize[2] = nz;

    //create the vectors for the representation of grid points
    for (unsigned int i = 0; i < nx; i++) {
        std::vector<std::vector<std::vector<unsigned int> > > currXVector;

        for (unsigned int j = 0; j < ny; j++) {
            std::vector<std::vector<unsigned int> > currYVector;

            for (unsigned int k = 0; k < nz; k++) {
                std::vector<unsigned int> currZVector;
                currYVector.push_back(currZVector);
            }
            currXVector.push_back(currYVector);
        }
        balloonList.push_back(currXVector);
    }
}

void StaticBalloonList::addCoordinateWithId(Coordinate3D input, unsigned int id) {

    if (objectAllocator.find(id) == objectAllocator.end()) {

        int u, v, w;
        Coordinate3D pos = input;

        u = (int) round((pos.x - lowerPoint.x) / gridConstant);
        v = (int) round((pos.y - lowerPoint.y) / gridConstant);
        w = (int) round((pos.z - lowerPoint.z) / gridConstant);

        if (u >= gridSize[0] || v >= gridSize[1] || w >= gridSize[2] ||
            u < 0 || v < 0 || w < 0) {
            ERROR_STDERR("sphere's position is out of balloonlist-boundary area. "
                         "Position:" <<
                                     pos.x);
            exit(1);
        } else {
            balloonList[u][v][w].push_back(id);
            coordinateAllocator[id] = input;
        }
    }

}

void StaticBalloonList::getInteractions(Coordinate3D myPos, std::vector<unsigned int> &neighbours) {

    std::vector<unsigned int>::iterator itIds;

    int u, v, w;
    Coordinate3D pos = myPos;

    u = round((pos.x - lowerPoint.x) / gridConstant);
    v = round((pos.y - lowerPoint.y) / gridConstant);
    w = round((pos.z - lowerPoint.z) / gridConstant);

    int nHSize;
    nHSize = (int) ceil(threshold / gridConstant);

    for (int i = u - nHSize; i <= u + nHSize; i++) {
        if (i < gridSize[0] && i >= 0) {
            for (int j = v - nHSize; j <= v + nHSize; j++) {
                if (j < gridSize[1] && j >= 0) {
                    for (int k = w - nHSize; k <= w + nHSize; k++) {
                        if (k < gridSize[2] && k >= 0) {
                            itIds = balloonList[i][j][k].begin();
                            while (itIds != balloonList[i][j][k].end()) {
                                unsigned int idPN = *itIds;
                                double distance = pos.calculateEuclidianDistance(coordinateAllocator[idPN]);
                                if (distance < threshold) {
                                    neighbours.push_back(idPN);
                                }
                                itIds++;
                            }
                        }
                    }
                }
            }
        }
    }

}

unsigned int StaticBalloonList::getClosestObjectIndex(Coordinate3D myPos) {
    std::vector<unsigned int>::iterator itIds;

    int u, v, w;
    Coordinate3D pos = myPos;

    u = round((pos.x - lowerPoint.x) / gridConstant);
    v = round((pos.y - lowerPoint.y) / gridConstant);
    w = round((pos.z - lowerPoint.z) / gridConstant);

    int nHSize;
    nHSize = (int) ceil(threshold / gridConstant);
    double minDistance = 1000;
    unsigned int closestObjectIndex = 999999;

    for (int i = u - nHSize; i <= u + nHSize; i++) {
        if (i < gridSize[0] && i >= 0) {
            for (int j = v - nHSize; j <= v + nHSize; j++) {
                if (j < gridSize[1] && j >= 0) {
                    for (int k = w - nHSize; k <= w + nHSize; k++) {
                        if (k < gridSize[2] && k >= 0) {
                            itIds = balloonList[i][j][k].begin();
                            while (itIds != balloonList[i][j][k].end()) {
                                unsigned int idPN = *itIds;
                                double distance = pos.calculateEuclidianDistance(coordinateAllocator[idPN]);
                                if (distance < minDistance) {
                                    minDistance = distance;
                                    closestObjectIndex = idPN;
                                }
                                itIds++;
                            }

                        }
                    }
                }
            }
        }
    }
    if (closestObjectIndex == 99999) {
        INFO_STDOUT("warning: no closest neighbour found!");
    }
    return closestObjectIndex;
}

void StaticBalloonList::getClosestObjectIndices(Coordinate3D myPos, std::vector<unsigned int> &neighbourList,
                                                unsigned int closestXParticles) {
    std::vector<unsigned int>::iterator itIds;
    std::vector<std::pair<double, unsigned int> > distancesOfIndices;

    int u, v, w;
    Coordinate3D pos = myPos;

    u = round((pos.x - lowerPoint.x) / gridConstant);
    v = round((pos.y - lowerPoint.y) / gridConstant);
    w = round((pos.z - lowerPoint.z) / gridConstant);

    int nHSize;
    nHSize = (int) ceil(threshold / gridConstant);
    double minDistance = 1000;
    unsigned int closestObjectIndex = 999999;

    for (int i = u - nHSize; i <= u + nHSize; i++) {
        if (i < gridSize[0] && i >= 0) {
            for (int j = v - nHSize; j <= v + nHSize; j++) {
                if (j < gridSize[1] && j >= 0) {
                    for (int k = w - nHSize; k <= w + nHSize; k++) {
                        if (k < gridSize[2] && k >= 0) {
                            itIds = balloonList[i][j][k].begin();
                            while (itIds != balloonList[i][j][k].end()) {
                                unsigned int idPN = *itIds;
                                double distance = pos.calculateEuclidianDistance(coordinateAllocator[idPN]);
                                std::pair<double, unsigned int> curDistanceOfIdx = std::pair<double, unsigned int>(
                                        distance, idPN);
                                distancesOfIndices.push_back(curDistanceOfIdx);
                                itIds++;
                            }
                        }
                    }
                }
            }
        }
    }

    std::sort(distancesOfIndices.begin(), distancesOfIndices.end());
    std::vector<std::pair<double, unsigned int> >::iterator it = distancesOfIndices.begin();
    unsigned int curNumberOfInserts = 0;
    while (it != distancesOfIndices.end() && curNumberOfInserts < closestXParticles) {
        neighbourList.push_back(it->second);
        curNumberOfInserts++;
        it++;
    }
}