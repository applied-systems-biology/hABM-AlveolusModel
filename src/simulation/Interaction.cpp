//  Copyright by Christoph Saffer, Paul Rudolph, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#include <sstream>

#include "simulation/Interaction.h"
#include "analyser/Analyser.h"
#include "simulation/InteractionState.h"
#include "simulation/Cell.h"
#include "simulation/InteractionStateFactory.h"
#include "simulation/InteractionFactory.h"


Interaction::Interaction(std::string identifier, Cell *cell1, Cell *cell2, double time_delta, double current_time)
        : interactionId(InteractionFactory::generateInteractionId()) {
    cellOne = cell1;
    cellTwo = cell2;
    identifier_ = identifier;
    isActiven = true;
    currentCondition = 0;
    setDelete = false;
}

void Interaction::setInitialState(double time_delta, double current_time, Cell *initiatingCell) {

    interactionState = InteractionStateFactory::createInteractionState(this, "InitialInteractionState", cellOne,cellTwo);
    if (cellularConditions.find(initiatingCell) != cellularConditions.end()) {
        currentCondition = cellularConditions[initiatingCell].get();
    }
    interactionState->stateTransition(time_delta, current_time);
}

void Interaction::setState(std::string nameOfState) {
    this->oldinteractionState = this->interactionState;
    this->interactionState = InteractionStateFactory::createInteractionState(this, nameOfState, cellOne, cellTwo);
}

void Interaction::handle(Cell *cell, double timestep, double current_time) {
    interactionState->handleInteraction(cell, timestep, current_time);
}

Cell *Interaction::getFirstCell() {
    return cellOne;
}

Cell *Interaction::getSecondCell() {
    return cellTwo;
}

Cell *Interaction::getOtherCell(Cell *cell) {
    Cell *otherCell;
    if (cell == cellOne) {
        otherCell = cellTwo;
    } else {
        otherCell = cellOne;
    }
    return otherCell;
}

bool Interaction::isActive() {
    return isActiven;
}

void Interaction::close() {
    isActiven = false;
}

void Interaction::fireInteractionEvent(InteractionEvent *ievent) {
    cellOne->recieveInteractionEvent(ievent);
    cellTwo->recieveInteractionEvent(ievent);
}

std::string Interaction::getInteractionName() const {
    return "Interaction";
}

void Interaction::includeInteractionXMLOutput(XMLFile *xmlFile, XMLNode *node, Cell *cell) {
    XMLNode interactionNode = xmlFile->addChildToNode(*node, this->getInteractionName());
    std::ostringstream ssid;
    ssid << this->getOtherCell(cell)->getId();
    xmlFile->addDataFieldToNode(interactionNode, "partner-id", "discrete", "unsigned int", ssid.str());
    xmlFile->addDataFieldToNode(interactionNode, "state", "discrete", "string", interactionState->getStateName());
}

Condition *Interaction::getCurrentCondition() {
    return currentCondition;
}

std::shared_ptr<Collision> Interaction::getNextCollision() {
    std::shared_ptr<Collision> coll = nullptr;
    if (!currentCollisions.empty()) {
        coll = currentCollisions.front();
        currentCollisions.pop();
    }
    return coll;
}