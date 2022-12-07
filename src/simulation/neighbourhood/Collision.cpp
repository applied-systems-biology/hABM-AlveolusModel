//  Copyright by Christoph Saffer, Paul Rudolph, Sandra Timme, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#include "Collision.h"
#include "simulation/Cell.h"
#include "simulation/AgentManager.h"


Collision::Collision(Cell *collisionCell, SphereRepresentation *mySphere, SphereRepresentation *collisionSphere,
                     double overlap) {
    this->cell = collisionCell->getSite()->getAgentManager()->getCellBySphereRepId(mySphere->getId());
    this->collisionCell = collisionCell;
    this->collisionSphere = collisionSphere;
    this->mySphere = mySphere;
    this->overlap = overlap;
}

void Collision::calculateTimeTillFirstContact() {
    double currentTimestep = cell->getMovement()->getCurrentTimestep();

    Coordinate3D v_a = Coordinate3D();
    v_a = *cell->getMovement()->getCurrentMove();
    v_a *= 1.0 / currentTimestep;

    double absV_a = v_a.getMagnitude();

    double radius_a, radius_1;
    radius_a = mySphere->getRadius();
    radius_1 = collisionSphere->getRadius();

    Coordinate3D r_1 = collisionSphere->getPosition();
    Coordinate3D r_a = cell->getPreviousPosition();
    Coordinate3D r_adt = mySphere->getPosition();

    Coordinate3D d_1a{r_a - r_1};
    double absD_1a = d_1a.getMagnitude();

    Coordinate3D termP = Coordinate3D();
    termP = v_a;
    double termP_result = termP.scalarProduct(d_1a);

    double p, q;
    p = (2.0 * termP_result) / (absV_a * absV_a);
    q = (absD_1a * absD_1a - (radius_1 + radius_a) * (radius_1 + radius_a)) / (absV_a * absV_a);

    double t1 = 0, t2 = 0, discr;
    discr = p * p / 4.0 - q;
    if (discr >= 0) {
        t1 = -p / 2.0 + sqrt(discr);
        t2 = -p / 2.0 - sqrt(discr);
    }

    if (t2 < 0) {
        timeToFirstContact = currentTimestep;
    } else {
        if (t2 < currentTimestep) {
            timeToFirstContact = t2;
        } else {
            timeToFirstContact = currentTimestep;
        }
    }

}