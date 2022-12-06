//  Copyright by Christoph Saffer, Paul Rudolph, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#ifndef AVOIDANCEINTERACTION_H
#define AVOIDANCEINTERACTION_H

#include <string>

#include "simulation/Interaction.h"

class AvoidanceInteraction : public Interaction {
public:
    // Class for interaction type by which cell are separated after interaction
    AvoidanceInteraction(std::string identifier, Cell *cell1, Cell *cell2, double time_delta, double current_time);
    [[nodiscard]] std::string getInteractionName() const final;
};

#endif    /* AVOIDANCEINTERACTION_H */

