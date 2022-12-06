//  Copyright by Christoph Saffer, Paul Rudolph, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#ifndef CELLSTATE_H
#define    CELLSTATE_H

#include <string>
#include <memory>
#include <vector>
#include <map>

#include "simulation/CellStateFactory.h"
#include "simulation/InteractionEvent.h"
#include "simulation/rates/ConstantRate.h"

class Cell;
class Randomizer;

class CellState {
public:
  // Class for wrapping cell states functionality
    CellState(const std::string &state_name, Cell *cell, std::map<std::string, const Rate *> next_states);

    ~CellState() = default;

    void handleInteractionEvent(InteractionEvent *interactionEvent);
    void stateTransition(double timestep, double current_time);
    std::string getStateName();
    bool checkForDeath(double current_time);
    void changeState(std::string stateName) {};
    void setNextState(std::string stateName);
    void addNextStateWithRate(const std::string &nameNextState, const Rate *rate);
    void addNextStateWithRate(const std::string &nameNextState, double rateOfNextState);

protected:
    void selectNextState(double timestep, Randomizer *randomizer);

    Cell *cell_;
    bool end_state_{};
    std::string next_state_;
    std::string current_state_;
    std::vector<std::unique_ptr<Rate>> own_rates_;
    std::map<std::string, const Rate *> next_states_rates_;
};

#endif    /* CELLSTATE_H */

