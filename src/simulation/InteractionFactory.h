//  Copyright by Christoph Saffer, Paul Rudolph, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#ifndef INTERACTIONFACTORY_H
#define    INTERACTIONFACTORY_H

#include <memory>
#include <utility>
#include <map>

#include "simulation/Interaction.h"
#include "analyser/Analyser.h"

class Collision;
class Cell;

class InteractionFactory {
public:
  // Factory class for the interactions specified in the simulerator configuration
    InteractionFactory() = delete;
    static void initialize(
            const std::vector<std::unique_ptr<abm::util::SimulationParameters::InteractionParameters>> &interaction_parameters);
    static void close();
    static unsigned int generateInteractionId();

    static std::shared_ptr<Interaction> createInteraction(double time_delta,
                                                          double current_time,
                                                          const std::shared_ptr<Collision> &collision,
                                                          InSituMeasurements *measurements);
    static std::shared_ptr<Interaction> createAvoidanceInteraction(Cell *cell_1,
                                                                   std::shared_ptr<Collision> collision,
                                                                   double time_delta,
                                                                   double current_time);
    static bool isInteractionsOn();

private:
    static std::tuple<std::string, std::string> retrieveInteractionIdentifier(Cell *cell_1, Cell *cell_2);
    static unsigned int interaction_id_;
    static std::map<std::string, std::string> interaction_types_;
    static std::map<std::string, std::vector<std::pair<std::string, std::vector<std::string>>>> interaction_conditions_;
    static std::map<std::pair<std::string, std::string>, std::string> interaction_pair_types_;

};

#endif    /* INTERACTIONFACTORY_H */

