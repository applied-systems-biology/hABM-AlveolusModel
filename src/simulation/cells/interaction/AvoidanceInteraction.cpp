//  Copyright by Christoph Saffer, Paul Rudolph, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#include "AvoidanceInteraction.h"
#include "simulation/InteractionState.h"


AvoidanceInteraction::AvoidanceInteraction(std::string identifier, Cell *cell1, Cell *cell2, double time_delta,
                                           double current_time) : Interaction(identifier, cell1, cell2, time_delta,
                                                                              current_time) {
    setInitialState(time_delta, current_time);
}

std::string AvoidanceInteraction::getInteractionName() const {
    return "AvoidanceInteraction";
}