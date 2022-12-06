//  Copyright by Christoph Saffer, Paul Rudolph, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#include "CellState.h"
#include "simulation/Cell.h"
#include "utils/macros.h"
#include "basic/Randomizer.h"


CellState::CellState(const std::string &state_name, Cell *cell, std::map<std::string, const Rate *> next_states) {
    current_state_ = state_name;
    next_states_rates_ = std::move(next_states);
    end_state_ = false;
    cell_ = cell;
}

std::string CellState::getStateName() {
    return current_state_;
}

void CellState::stateTransition(double timestep, double current_time) {
    selectNextState(timestep, cell_->getSite()->getRandomGenerator());

    if (next_state_ == "self") {
        //in principle do nothing
    } else {
        //change the state and fire event
        auto nextCellState = cell_->getCellStateByName(next_state_);
        if (nextCellState != 0) {
            cell_->setState(nextCellState);
        } else {
            ERROR_STDERR("there is a problem with CellState changes->" +
                         getStateName() + " ns->" + next_state_);
        }
        cell_->setTimestepLastTreatment(current_time);
    }

    if (end_state_) {
        checkForDeath(current_time);
    }
    next_state_ = "";
}

void CellState::handleInteractionEvent(InteractionEvent *interactionEvent) {

}

bool CellState::checkForDeath(double current_time) {
    return current_state_ == "Death";
}

void CellState::selectNextState(double timestep, Randomizer *randomizer) {

    if (next_states_rates_.empty()) {
        next_state_ = "self";
    } else {
        double bottom = 0;
        double top = 0;
        std::string backup{};
        double p = randomizer->generateDouble();
        for (const auto&[kRateName, kCurRate]: next_states_rates_) {
            double cur_prob = kCurRate->calculateProbability(timestep, nullptr, cell_, this->cell_->getSite());
            if (cur_prob < 0) {
                backup = kRateName;
            } else {
                top += cur_prob;
                if (p >= bottom && p < top) {
                    next_state_ = kRateName;
                    break;
                    //next state was successfully found
                }
                bottom = top;
            }
        }
        if (next_state_.empty()) {
            if (backup.empty()) {
                next_state_ = "self";
            } else {
                next_state_ = backup;
            }
        }
    }
}

void CellState::addNextStateWithRate(const std::string &nameNextState, const double rateOfNextState) {
    const auto &rate = own_rates_.emplace_back(std::make_unique<ConstantRate>(rateOfNextState));
    next_states_rates_[nameNextState] = rate.get();
}

void CellState::addNextStateWithRate(const std::string &nameNextState, const Rate *rate) {
    next_states_rates_[nameNextState] = rate;
}

void CellState::setNextState(std::string stateName) {
    next_state_ = stateName;
}
