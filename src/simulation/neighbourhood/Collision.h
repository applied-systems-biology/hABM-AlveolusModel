//  Copyright by Christoph Saffer, Paul Rudolph, Sandra Timme, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#ifndef COLLISION_H
#define    COLLISION_H

#include <map>
#include <set>

#include "simulation/morphology/SphereRepresentation.h"

class Cell;

enum class MeasurementType {
    NO_INTERACTION,
    EXISTING_INTERACTION,
    NEW_INTERACTION
};

class Collision {
public:
  // Class for collision check between cells represented as spheres
    Collision() = default;
    Collision(Cell *collisionCell, SphereRepresentation *mySphere, SphereRepresentation *collisionSphere,
              double overlap);
    [[nodiscard]] Cell *getCell() const { return cell; };
    [[nodiscard]] Cell *getCollisionCell() const { return collisionCell; };
    SphereRepresentation *getMySphere() { return mySphere; };
    SphereRepresentation *getCollisionSphere() { return collisionSphere; };
    double getTimeToFirstContact() { return timeToFirstContact; };
    void calculateTimeTillFirstContact();

    MeasurementType type{MeasurementType::NO_INTERACTION};
private:
    Cell *cell{};
    Cell *collisionCell{};
    SphereRepresentation *mySphere{};
    SphereRepresentation *collisionSphere{};
    double timeToFirstContact{};
    double overlap{};
};

#endif    /* COLLISION_H */

