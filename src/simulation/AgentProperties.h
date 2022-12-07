//  Copyright by Christoph Saffer, Paul Rudolph, Sandra Timme, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#ifndef AGENTPROPERTIES_H
#define    AGENTPROPERTIES_H

#include <iostream>

#include "simulation/Movement.h"
#include "simulation/Morphology.h"

class Interactions;
class AgentProperties {

public:
  // Class for separating some default properties of agents from the agent class.
    AgentProperties() = default;

    void setMorphology(std::shared_ptr<Morphology> morphology);
    void setMovement(std::shared_ptr<Movement> movement);
    void setPassiveMovement(std::shared_ptr<Movement> passiveMovement);
    void setInteractions(std::shared_ptr<Interactions> interactions);

    Movement *getMovement();
    Movement *getPassiveMovement();
    Interactions *getInteractions();
    Morphology *getMorphology();
private:
    std::shared_ptr<Movement> movement_{};
    std::shared_ptr<Movement> passive_movement_{};
    std::shared_ptr<Interactions> interactions_{};
    std::shared_ptr<Morphology> morphology_{};
};

#endif    /* AGENTPROPERTIES_H */

