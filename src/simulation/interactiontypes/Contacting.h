//  Copyright by Christoph Saffer, Paul Rudolph, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#ifndef CONTACTING_H
#define    CONTACTING_H

#include "simulation/InteractionType.h"

class Contacting : public InteractionType {
public:
  // Class for 'Contacting' interaction that is the default interaction between similar cells on collision
    Contacting(bool adhere, double mustOverhead) : adhere(adhere), mustOverhead(mustOverhead) {}

    void handleInteraction(Interaction *interaction, Cell *cell, double timestep, double current_time) final;
    [[nodiscard]] std::string getTypeName() const final;

private:
    bool adhere;
    double mustOverhead;
};

#endif    /* CONTACTING_H */

