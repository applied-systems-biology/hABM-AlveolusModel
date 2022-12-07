//  Copyright by Christoph Saffer, Paul Rudolph, Sandra Timme, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#include "simulation/AgentProperties.h"
#include "simulation/Interactions.h"

void AgentProperties::setMovement(std::shared_ptr<Movement> movement) {
    movement_ = std::move(movement);
}

void AgentProperties::setPassiveMovement(std::shared_ptr<Movement> passiveMovement) {
    passive_movement_ = std::move(passiveMovement);
}

Movement *AgentProperties::getMovement() {
    return movement_.get();
}

Movement *AgentProperties::getPassiveMovement() {
    return passive_movement_.get();
}

void AgentProperties::setInteractions(std::shared_ptr<Interactions> interactions) {
    interactions_ = interactions;
}

Interactions *AgentProperties::getInteractions() {
    return interactions_.get();
}

void AgentProperties::setMorphology(std::shared_ptr<Morphology> morphology) {
    morphology_ = std::move(morphology);
}

Morphology *AgentProperties::getMorphology() {
    if (morphology_ == nullptr) {
        morphology_ = std::make_unique<Morphology>();
    }

    return morphology_.get();
}