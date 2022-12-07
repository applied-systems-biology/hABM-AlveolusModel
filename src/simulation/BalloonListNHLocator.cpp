//  Copyright by Christoph Saffer, Paul Rudolph, Sandra Timme, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#include <cmath>
#include <cstdlib>

#include "simulation/neighbourhood/Collision.h"
#include "simulation/BalloonListNHLocator.h"
#include "utils/macros.h"
#include "simulation/AgentManager.h"
#include "simulation/Site.h"

BalloonListNHLocator::BalloonListNHLocator(double gridConstant, Coordinate3D lowerValues, Coordinate3D upperValues,
                                           Site *site) : NeighbourhoodLocator(site) {
    this->gridConstant = gridConstant;
    this->lowerPoint = lowerValues;
    this->upperPoint = upperValues;
    complexes = 0;
    checksum++;
    instantiate(site);
}

BalloonListNHLocator::~BalloonListNHLocator() {
}

void BalloonListNHLocator::instantiate(Site *site) {
    initialGridCreation();
}

std::vector<std::shared_ptr<Collision> > BalloonListNHLocator::getCollisions(Agent *agent) {

    std::vector<SphereRepresentation *>::iterator itCellsSpheres;
    std::vector<SphereRepresentation *> currentCellsSpheres;
    currentCellsSpheres = agent->getAgentProperties()->getMorphology()->getAllSpheresOfThis();
    itCellsSpheres = currentCellsSpheres.begin();
    std::vector<std::shared_ptr<Collision>> neighbours;
    SphereRepresentation *currentCellsSphere;

    while (itCellsSpheres != currentCellsSpheres.end()) {

        currentCellsSphere = *itCellsSpheres;
        std::vector<int> agentGridPoint(3);
        int u, v, w;

        agentGridPoint = sphereRepresentationAllocator[currentCellsSphere];

        u = agentGridPoint[0];
        v = agentGridPoint[1];
        w = agentGridPoint[2];

        int nHSize;
        nHSize = ceil(thresholdDistance / gridConstant);

        for (int i = u - nHSize; i <= u + nHSize; i++) {
            if (i >= 0 && i < (gridSize[0])) {
                for (int j = v - nHSize; j <= v + nHSize; j++) {
                    if (j >= 0 && j < (gridSize[1])) {
                        for (int k = w - nHSize; k <= w + nHSize; k++) {
                            if (k >= 0 && k < (gridSize[2])) {
                                checkCollisions(&neighbours, agent, currentCellsSphere, i, j, k);
                            }
                        }
                    }
                }
            }
        }
        itCellsSpheres++;
    }
    return neighbours;
}

void BalloonListNHLocator::updateDataStructures(SphereRepresentation *sphereRep) {

    if (sphereRepresentationAllocator.find(sphereRep) != sphereRepresentationAllocator.end()) {
        std::vector<int> sphereGridPoint(3);
        sphereGridPoint = sphereRepresentationAllocator[sphereRep];
        unsigned int uOld, vOld, wOld;
        uOld = sphereGridPoint[0];
        vOld = sphereGridPoint[1];
        wOld = sphereGridPoint[2];

        unsigned int u, v, w;
        Coordinate3D pos = sphereRep->getPosition();

        u = (unsigned int) round((pos.x - lowerPoint.x) / gridConstant);
        v = (unsigned int) round((pos.y - lowerPoint.y) / gridConstant);
        w = (unsigned int) round((pos.z - lowerPoint.z) / gridConstant);

        if (uOld == u && vOld == v && wOld == w) {
            //do nothing
        } else {
            removeSphereRepresentation(sphereRep);
            addSphereRepresentation(sphereRep);
        }
    } else {
        DEBUG_STDOUT("The sphere is not yet in the system, can not do update");
    }
}

void BalloonListNHLocator::removeSphereRepresentation(SphereRepresentation *sphereRep) {
    std::vector<SphereRepresentation *>::iterator toDelete;

    std::vector<int> sphereGridPoint(3);
    if (sphereRepresentationAllocator.find(sphereRep) != sphereRepresentationAllocator.end()) {
        sphereGridPoint = sphereRepresentationAllocator[sphereRep];
        unsigned int u, v, w;
        u = sphereGridPoint[0];
        v = sphereGridPoint[1];
        w = sphereGridPoint[2];

        sphereRepresentationAllocator.erase(sphereRep);

        toDelete = remove(balloonList[u][v][w].begin(), balloonList[u][v][w].end(), sphereRep);
        balloonList[u][v][w].erase(toDelete, balloonList[u][v][w].end());
    }

}

void BalloonListNHLocator::addSphereRepresentation(SphereRepresentation *sphereRep) {
    if (sphereRepresentationAllocator.find(sphereRep) == sphereRepresentationAllocator.end()) {
        std::vector<int> position(3);

        int u, v, w;
        Coordinate3D pos = sphereRep->getPosition();

        u = (int) round((pos.x - lowerPoint.x) / gridConstant);
        v = (int) round((pos.y - lowerPoint.y) / gridConstant);
        w = (int) round((pos.z - lowerPoint.z) / gridConstant);

        position[0] = u;
        position[1] = v;
        position[2] = w;

        if (site_->getNumberOfSpatialDimensions() == 2) {
            if (u >= gridSize[0] || v >= gridSize[1] || u < 0 || v < 0) {
                ERROR_STDERR("Sphere's position is out of balloonlist-boundary area. "
                             "Position: " << pos.printCoordinates());
                exit(1);
            } else {
                balloonList[u][v][w].push_back(sphereRep);

                sphereRepresentationAllocator[sphereRep] = position;
            }
        } else if (site_->getNumberOfSpatialDimensions() == 3) {
            if (u >= gridSize[0] || v >= gridSize[1] || w >= gridSize[2] ||
                u < 0 || v < 0 || w < 0) {
                ERROR_STDERR("Sphere's position is out of balloonlist-boundary area. "
                             "Position: (" << pos.x << ", " << pos.y << ", " << pos.z << ")");

                exit(1);
            } else {
                balloonList[u][v][w].push_back(sphereRep);

                sphereRepresentationAllocator[sphereRep] = position;
            }
        }
    }
}

void BalloonListNHLocator::initialGridCreation() {
    int nx, ny, nz;
    nx = (int) ceil((upperPoint.x - lowerPoint.x) / gridConstant) + 1;
    ny = (int) ceil((upperPoint.y - lowerPoint.y) / gridConstant) + 1;
    nz = (int) ceil((upperPoint.z - lowerPoint.z) / gridConstant) + 1;

    gridSize[0] = nx;
    gridSize[1] = ny;
    gridSize[2] = nz;

    //create the vectors for the representation of grid points
    for (int i = 0; i < nx; i++) {
        std::vector<std::vector<std::vector<SphereRepresentation *> > > currXVector;
        for (int j = 0; j < ny; j++) {
            std::vector<std::vector<SphereRepresentation *> > currYVector;
            for (int k = 0; k < nz; k++) {
                std::vector<SphereRepresentation *> currZVector;
                currZVector.reserve(10);
                currYVector.push_back(currZVector);
            }
            currXVector.push_back(currYVector);
        }
        balloonList.push_back(currXVector);
    }
}

bool BalloonListNHLocator::checkCollisions(std::vector<std::shared_ptr<Collision> > *collisions, Agent *agent,
                                           SphereRepresentation *sphereRep, int u, int v, int w, bool justCheck) {
    bool returnVal = false;

    double distance, minDistance, r1, r2;
    for (auto currNeighbour: balloonList[u][v][w]) {
        Cell *collisionCell = agent->getSite()->getAgentManager()->getCellBySphereRepId(currNeighbour->getId());
        if (collisionCell != nullptr) {
            if (collisionCell->getId() != agent->getId() && !collisionCell->isDeleted()) {

                distance = sphereRep->getPosition().calculateEuclidianDistance(currNeighbour->getPosition());
                r1 = sphereRep->getRadius();
                r2 = currNeighbour->getRadius();
                minDistance = r1 + r2;

                if (distance <= minDistance) {
                    returnVal = true;
                    if (!justCheck) {
                        collisions->push_back(std::make_shared<Collision>(collisionCell, sphereRep, currNeighbour,
                                                            (minDistance - distance)));
                    }
                }
            }
        }
    }
    return returnVal;
}

bool BalloonListNHLocator::hasCollision(Agent *agent) {
    bool hasOneCollision = false;
    std::string agentName = agent->getTypeName();
    std::vector<SphereRepresentation *>::iterator itCellsSpheres;
    std::vector<SphereRepresentation *> currentCellsSpheres;

    currentCellsSpheres = agent->getAgentProperties()->getMorphology()->getAllSpheresOfThis();

    itCellsSpheres = currentCellsSpheres.begin();
    std::vector<std::shared_ptr<Collision>> neighbours;
    SphereRepresentation *currentCellsSphere;

    while (itCellsSpheres != currentCellsSpheres.end()) {

        currentCellsSphere = *itCellsSpheres;
        std::vector<int> agentGridPoint(3);
        int u, v, w;

        agentGridPoint = sphereRepresentationAllocator[currentCellsSphere];

        u = agentGridPoint[0];
        v = agentGridPoint[1];
        w = agentGridPoint[2];

        int nHSize;
        nHSize = ceil(thresholdDistance / gridConstant);

        for (int i = u - nHSize; i <= u + nHSize; i++) {
            if (i >= 0 && i < ((int) gridSize[0])) {
                for (int j = v - nHSize; j <= v + nHSize; j++) {
                    if (j >= 0 && j < ((int) gridSize[1])) {
                        for (int k = w - nHSize; k <= w + nHSize; k++) {
                            if (k >= 0 && k < ((int) gridSize[2])) {
                                if ((hasOneCollision = checkCollisions(&neighbours, agent, currentCellsSphere, i, j, k, true))) {
                                    break;
                                }
                            }
                        }
                    }
                    if (hasOneCollision) {
                        break;
                    }
                }
            }
            if (hasOneCollision) {
                break;
            }
        }
        if (hasOneCollision) {
            break;
        }
        itCellsSpheres++;
    }

    return hasOneCollision;
}

std::vector<Coordinate3D>
BalloonListNHLocator::getCollisionSpheres(SphereRepresentation *sphereRep, Coordinate3D dirVec) {
    std::vector<Coordinate3D> collisionPositions;
    std::vector<std::shared_ptr<Collision>> neighbours;
    SphereRepresentation *currentCellsSphere = sphereRep;
    std::vector<int> agentGridPoint(3);

    double newl = ((currentCellsSphere->getRadius()) / dirVec.getMagnitude());
    Coordinate3D sphpos = currentCellsSphere->getPosition();
    Coordinate3D futpos = {sphpos.x + newl * dirVec.x, sphpos.y + newl * dirVec.y, sphpos.z + newl * dirVec.z};

    int u, v, w;

    agentGridPoint = sphereRepresentationAllocator[currentCellsSphere];

    u = (int) agentGridPoint[0];
    v = (int) agentGridPoint[1];
    w = (int) agentGridPoint[2];

    int nHSize;
    nHSize = (int) ceil(thresholdDistance / gridConstant);

    for (int i = u - nHSize; i <= u + nHSize; i++) {
        if (i >= 0 && i < ((int) gridSize[0])) {
            for (int j = v - nHSize; j <= v + nHSize; j++) {
                if (j >= 0 && j < ((int) gridSize[1])) {
                    for (int k = w - nHSize; k <= w + nHSize; k++) {
                        if (k >= 0 && k < ((int) gridSize[2])) {
                            for (auto x: balloonList[i][j][k]) {
                                double min_dist = (x->getRadius() + currentCellsSphere->getRadius());
                                double dist = x->getPosition().calculateEuclidianDistance(futpos);
                                if (min_dist > dist) {
                                    collisionPositions.emplace_back(x->getPosition());
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return collisionPositions;
}

int BalloonListNHLocator::controlFunction() {
    int count = 0;
    for (int i = 0; i < ((int) gridSize[0]); i++) {
        for (int j = 0; j < ((int) gridSize[1]); j++) {
            for (int k = 0; k < ((int) gridSize[2]); k++) {
                count += balloonList[i][j][k].size();
            }
        }
    }
    return count;
}

std::string BalloonListNHLocator::getTypeName() {
    return "BalloonListNHLocator";
}

int BalloonListNHLocator::getNumberOfAgentTypeInBalloonList(std::string agentType) {
    int count = 0;
    for (size_t i1 = 0; i1 < balloonList.size(); i1++) {
        for (size_t i2 = 0; i2 < balloonList[i1].size(); i2++) {
            for (size_t i3 = 0; i3 < balloonList[i1][i2].size(); i3++) {
                for (size_t i4 = 0; i4 < balloonList[i1][i2][i3].size(); i4++) {
                    std::string type = balloonList[i1][i2][i3][i4]->getMorphologyElementThisBelongsTo()->getMorphologyThisBelongsTo()->getCellThisBelongsTo()->getTypeName();
                    if (type == agentType) {
                        count++;
                    }
                }
            }
        }
    }
    return count;
}