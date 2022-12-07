
//  Copyright by Christoph Saffer, Paul Rudolph, Sandra Timme, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#include "simulation/InteractionEvent.h"


InteractionEvent::InteractionEvent(std::string previousState, std::string nextState) {
    this->previousState = previousState;
    this->nextState = nextState;
    setDescriptiveName();
}

InteractionEvent::InteractionEvent(std::string previousState, std::string nextState, Interaction *interaction) {
    this->previousState = previousState;
    this->nextState = nextState;
    this->interaction = interaction;
    setDescriptiveName();
}

void InteractionEvent::setDescriptiveName() {

}

std::string InteractionEvent::getDescriptiveName() {
    return descriptiveName;
}

std::string InteractionEvent::getNextState() {
    return nextState;
}

std::string InteractionEvent::getPreviousState() {
    return previousState;
}

Interaction *InteractionEvent::getInteraction() {
    return interaction;
}