//  Copyright by Christoph Saffer, Paul Rudolph, Sandra Timme, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#include "simulation/InteractionStateFactory.h"
#include "PhagocyteFungusInteraction.h"
#include "analyser/Analyser.h"
#include "simulation/Cell.h"


PhagocyteFungusInteraction::PhagocyteFungusInteraction(std::string identifier,
                                                       Cell *cell1,
                                                       Cell *cell2,
                                                       double time_delta,
                                                       double current_time) : Interaction(identifier,
                                                                                          cell1,
                                                                                          cell2,
                                                                                          time_delta,
                                                                                          current_time) {
    auto phagoCond = std::make_shared<Condition>(cellOne);
    auto fungusCond = std::make_shared<Condition>(cellTwo);
    cellularConditions[cellOne] = phagoCond;
    cellularConditions[cellTwo] = fungusCond;
    cellularStringConditions[cellOne->getTypeName()] = phagoCond;
    cellularStringConditions[cellTwo->getTypeName()] = fungusCond;
    interactionState = InteractionStateFactory::createInteractionState(this, "InitialInteractionState", cellOne,
                                                                       cellTwo);
    if (cellularConditions.find(cell1) != cellularConditions.end()) {
        currentCondition = cellularConditions[cell1].get();
    }
}

PhagocyteFungusInteraction::PhagocyteFungusInteraction(std::string identifier,
                                                       Cell *cell1,
                                                       Cell *cell2,
                                                       bool noInitialSetup,
                                                       double time_delta,
                                                       double current_time) : Interaction(identifier,
                                                                                          cell1,
                                                                                          cell2,
                                                                                          time_delta,
                                                                                          current_time) {
    auto phagoCond = std::make_shared<Condition>(cellOne);
    auto fungusCond = std::make_shared<Condition>(cellTwo);
    cellularConditions[cellOne] = phagoCond;
    cellularConditions[cellTwo] = fungusCond;
    cellularStringConditions[cellOne->getTypeName()] = phagoCond;
    cellularStringConditions[cellTwo->getTypeName()] = fungusCond;

    if (cellularConditions.find(cell1) != cellularConditions.end()) {
        currentCondition = cellularConditions[cell1].get();
    }
}

void PhagocyteFungusInteraction::handle(Cell *cell, double timestep, double current_time) {
    currentCondition = cellularConditions[cell].get();
    this->Interaction::handle(cell, timestep, current_time);
}

std::string PhagocyteFungusInteraction::getInteractionName() const {
    return "PhagocyteFungusInteraction";
}