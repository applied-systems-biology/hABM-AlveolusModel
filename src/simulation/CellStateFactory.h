//  Copyright by Christoph Saffer, Paul Rudolph, Sandra Timme, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#ifndef CELLSTATEFACTORY_H
#define CELLSTATEFACTORY_H

#include "utils/io_util.h"
#include "external/xmlParser/xmlParser.h"

class CellState;
class Cell;
class Rate;

class CellStateFactory {
    using StateSetup = std::map<std::string, std::map<std::string, const Rate *>>;
public:
  // Factory class for initializing all possible cell states according to the simulator configuration
    CellStateFactory() = delete;

    static void initialize(const std::unique_ptr<abm::util::SimulationParameters::SiteParameters> &site_parameters);
    static void close();
    static std::shared_ptr<CellState> createCellState(Cell *cell, const std::string &state_name);
    static const StateSetup &getStateSetup(const std::string &site_identifier, const std::string &agent_type);

private:
    static std::map<std::string, const Rate *> getNextStates(const std::string &stae_name,
                                                             const std::string &site_identifier,
                                                             const std::string &agent_type);
    static std::map<std::pair<std::string, std::string>, StateSetup> next_states;

};

#endif    /* CELLSTATEFACTORY_H */

