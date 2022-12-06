//  Copyright by Christoph Saffer, Paul Rudolph, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#ifndef RIGIDCONTACTING_H
#define    RIGIDCONTACTING_H

#include "simulation/InteractionType.h"
#include "simulation/Cell.h"

class RigidContacting : public InteractionType {
public:
    RigidContacting(double mustOverhead) : mustOverhead(mustOverhead) {}
    void handleInteraction(Interaction *interaction, Cell *cell, double timestep, double current_time) final;
    std::string getTypeName() const final { return "RigidContacting"; };

private:
    double mustOverhead;
};

#endif    /* RIGIDCONTACTING_H */

