//  Copyright by Christoph Saffer, Paul Rudolph, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#ifndef INTERACTIONTYPE_H
#define    INTERACTIONTYPE_H

#include <string>

class Cell;
class Interaction;

class InteractionType {
public:
  //  Default class for describing the interaction type of a cell-cell interaction.
    InteractionType() = default;
    virtual ~InteractionType() = default;
    virtual void handleInteraction(Interaction *interaction, Cell *cell, double timestep, double current_time) {};
    [[nodiscard]] virtual std::string getTypeName() const { return "InteractionType"; };
protected:
};

#endif    /* INTERACTIONTYPE_H */

