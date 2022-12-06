//  Copyright by Christoph Saffer, Paul Rudolph, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#include "simulation/InteractionState.h"

#include "simulation/Interaction.h"
#include "simulation/InteractionEvent.h"
#include "analyser/InSituMeasurements.h"
#include "simulation/Site.h"


void InteractionState::handleInteraction(Cell *cell, double timestep, double current_time) {
    interaction_type_->handleInteraction(interaction_, cell, timestep, current_time);
    if (!interaction_->isDelted()) {
        stateTransition(timestep, current_time);
    }
}

void InteractionState::stateTransition(double timestep, double current_time) {
    Cell *cell1 = interaction_->getFirstCell();
    Cell *cell2 = interaction_->getSecondCell();
    Condition *curCondition = interaction_->getCurrentCondition();
    selectNextState(timestep, cell1->getSite()->getRandomGenerator(), curCondition);

    if (cell1->agentTreatedInCurrentTimestep(current_time) || cell2->agentTreatedInCurrentTimestep(current_time)) {
        if (next_state_ != "NoInterplay" && next_state_.compare("self") != 0 && next_state_.compare("Avoidance") != 0) {
            auto iterNextStates = next_states_rates_.begin();
            while (iterNextStates != next_states_rates_.end()) {
                std::string nextStateName = iterNextStates->first;
                if (nextStateName == "NoInterplay") {
                    next_state_ = "NoInterplay";
                    break;
                }
                iterNextStates++;
            }
            if (next_state_ != "NoInterplay") {
                next_state_ = "self";
            }
        }
    }

    if (next_state_.compare("self") != 0 && next_state_ != "NoInterplay" && next_state_.compare("Avoidance") != 0) {
        cell1->setTimestepLastTreatment(current_time);
        cell2->setTimestepLastTreatment(current_time);
    }

    if (next_state_ != "self") {
        interaction_->setState(next_state_);
        if (interaction_->getCurrentState()->getStateName() == "Phagocytose") {
            //                cout << "[InteractionState] next state is
            //                Phagocytose" << '\n';
            Cell *cell1 = interaction_->getFirstCell();
            Cell *cell2 = interaction_->getSecondCell();
            cell1->getSite()->getMeasurments()->increment<PairMeasurement>("Phagocytosis-NC", "NCPhag", 1);
            cell1->getSite()->getMeasurments()->increment<PairMeasurement>("Phagocytosis-MC", "MCPhag", 1);
        }
        const auto ievent = std::make_unique<InteractionEvent>(getStateName(), next_state_, interaction_);
        interaction_->fireInteractionEvent(ievent.get());
    }

    if (interaction_->getCurrentState()->isEndState()) {
        interaction_->close();
    }
    next_state_ = "";
}

void InteractionState::fireInteractionEvent(const std::string &nextState) {
    const auto ievent = std::make_unique<InteractionEvent>(getStateName(), nextState, interaction_);
    interaction_->fireInteractionEvent(ievent.get());
}

void InteractionState::addNextStateWithRate(const std::string &nameNextState, const Rate *rate) {
    next_states_rates_[nameNextState] = rate;
}

void InteractionState::addNextStateWithRate(std::map<std::string, const Rate *> next_states_rates) {
    next_states_rates_ = next_states_rates;
}

bool InteractionState::isEndState() const { return end_state_; }

void InteractionState::selectNextState(double timestep, Randomizer *randomizer, Condition *condition) {
    if (next_states_rates_.empty()) {
        next_state_ = "self";
    } else {
        double bottom = 0;
        double top = 0;
        std::string backup{};
        double p = randomizer->generateDouble();
        for (const auto&[kRateName, kCurRate]: next_states_rates_) {
            double cur_prob = kCurRate->calculateProbability(timestep, condition, nullptr,nullptr);
            if (cur_prob < 0) {
                backup = kRateName;
            } else {
                top += cur_prob;
                if (p >= bottom && p < top) {
                    next_state_ = kRateName;
                    break;
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

const std::string &InteractionState::getStateName() const {
    return current_state_;
}

std::string InteractionState::getInteractionType() const {
    return interaction_type_->getTypeName();
}
