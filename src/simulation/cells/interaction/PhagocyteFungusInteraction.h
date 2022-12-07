//  Copyright by Christoph Saffer, Paul Rudolph, Sandra Timme, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#ifndef PHAGOCYTEFUNGUSINTERACTION_H
#define    PHAGOCYTEFUNGUSINTERACTION_H

#include <string>

#include "simulation/Interaction.h"


class PhagocyteFungusInteraction : public Interaction {

public:
  // Class for phagocytosis interaction. This class provides the main functionality if a phagocytosis event is triggered in a event chain.
    PhagocyteFungusInteraction(std::string identifier,
                               Cell *cellOne,
                               Cell *cellTwo,
                               double time_delta,
                               double current_time);
    PhagocyteFungusInteraction(std::string identifier,
                               Cell *cell1,
                               Cell *cell2,
                               bool noInitialSetup,
                               double time_delta,
                               double current_time);

    void handle(Cell *cell, double timestep, double current_time) final;
    [[nodiscard]] std::string getInteractionName() const final;
};

#endif    /* PHAGOCYTEFUNGUSINTERACTION_H */

