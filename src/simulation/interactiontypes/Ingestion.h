//  Copyright by Christoph Saffer, Paul Rudolph, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#ifndef INGESTION_H
#define    INGESTION_H

#include "simulation/InteractionType.h"


class Ingestion : public InteractionType {
public:
  // Class for interaction type 'Ingestion'. This class provides functionality for phagocytosis.
    Ingestion() = default;
    void handleInteraction(Interaction *interaction, Cell *cell, double timestep, double current_time) final;
    std::string getTypeName() const final;

private:
    bool isIngestingType(Cell *cell);
};

#endif    /* INGESTION_H */

