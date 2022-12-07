//  Copyright by Christoph Saffer, Paul Rudolph, Sandra Timme, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#ifndef ASPERGILLUSFUMIGATUS_H
#define    ASPERGILLUSFUMIGATUS_H

#include "simulation/Cell.h"

class AspergillusFumigatus : public Cell {
public:
  // Class for fungal cell that specify properties such as the information that a cell is phagocytosed.
    AspergillusFumigatus(std::unique_ptr<Coordinate3D> c, int id, Site *site, double time_delta, double current_time)
            : Cell(std::move(c),
                   id,
                   site,
                   time_delta,
                   current_time) {}

    void move(double timestep, double current_time) final;
    void doMorphologicalChanges(double timestep, double current_time) final;
    std::string getTypeName() final;
    void setup(double time_delta, double current_time, abm::util::SimulationParameters::AgentParameters *parameters) final;

private:
    void handleInteractionEvent(InteractionEvent *ievent) final;
    bool phagocytosed{};
};

#endif    /* ASPERGILLUSFUMIGATUS_H */

