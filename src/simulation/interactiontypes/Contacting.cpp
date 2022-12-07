//  Copyright by Christoph Saffer, Paul Rudolph, Sandra Timme, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#include "Contacting.h"
#include "simulation/Interaction.h"
#include "simulation/Cell.h"


std::string Contacting::getTypeName() const{
    return "Contacting";
}

void Contacting::handleInteraction(Interaction *interaction,Cell* cell, double timestep, double current_time){
    Cell *passiveCell,*activeCell;
    SphereRepresentation *passiveSphere, *activeSphere;
    
    std::shared_ptr<Collision> currentCollision;
    if (adhere) {
        passiveCell = cell;
        activeCell = interaction->getOtherCell(cell);
        activeSphere = activeCell->getAgentProperties()->getMorphology()->getBasicSphereOfThis();
        passiveSphere = passiveCell->getAgentProperties()->getMorphology()->getBasicSphereOfThis();
        mustOverhead = 0.0;
    } else {
        passiveCell = interaction->getOtherCell(cell);
        activeCell = cell;
    }
    int collisionNo=0;
    while ((currentCollision = interaction->getNextCollision()) != 0) {
        collisionNo++;
        if (adhere) {
            passiveSphere = currentCollision->getMySphere();
            activeSphere = currentCollision->getCollisionSphere();
            mustOverhead = 0.0;
        } else {
            passiveSphere = passiveCell->getAgentProperties()->getMorphology()->getBasicSphereOfThis();
            activeSphere = currentCollision->getMySphere();;
        }

        Coordinate3D backShift = cell->getSite()->generateBackShiftOnContacting(activeSphere,passiveSphere,mustOverhead);
        if(!cell->isDeleted()){
            activeCell->shiftPosition(&backShift,current_time, activeSphere,getTypeName());
        }

        currentCollision = nullptr;
    }   
}