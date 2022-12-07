//  Copyright by Christoph Saffer, Paul Rudolph, Sandra Timme, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#include <variant>
#include <utils/misc_util.h>

#include "simulation/InteractionStateFactory.h"
#include "io/InputConfiguration.h"
#include "simulation/interactiontypes/Ingestion.h"
#include "simulation/interactiontypes/Contacting.h"
#include "simulation/interactiontypes/RigidContacting.h"
#include "simulation/Interaction.h"
#include "simulation/Rate.h"
#include "utils/macros.h"
#include "simulation/RateFactory.h"

std::map<std::string, InteractionStateFactory::StateSetup>InteractionStateFactory::state_parameters_{};

void InteractionStateFactory::initialize(
        const std::vector<std::unique_ptr<abm::util::SimulationParameters::InteractionParameters>> &interaction_parameters) {
    if (RateFactory::isInitialized()) {
        for (const auto &parameters: interaction_parameters) {
            for (const auto &state: parameters->states) {
                std::map<std::string, const Rate *> state_setup;
                for (const auto&[next_state, rate_name]:state.next_states) {
                    state_setup.emplace(next_state, RateFactory::getRate(rate_name));
                }
                if (state.interaction_type == "InteractionType") {
                    state_parameters_[parameters->name].emplace(state.name, std::make_pair(InteractionType(),
                                                                                           std::move(state_setup)));
                } else if (state.interaction_type == "Contacting") {
                    state_parameters_[parameters->name].emplace(state.name,
                                                                std::make_pair(
                                                                        Contacting(state.adhere, state.must_overhead),
                                                                        std::move(state_setup)));
                } else if (state.interaction_type == "RigidContacting") {
                    state_parameters_[parameters->name].emplace(state.name,
                                                                std::make_pair(RigidContacting(state.must_overhead),
                                                                               std::move(state_setup)));
                } else if (state.interaction_type == "Ingestion") {
                    state_parameters_[parameters->name].emplace(state.name,
                                                                std::make_pair(Ingestion(), std::move(state_setup)));
                }
            }
        }
    } else {
        ERROR_STDERR("Rate Factory needs to be initialized before InteractionStateFactory.");
        exit(1);
    }
}

void InteractionStateFactory::close() {
    state_parameters_.clear();
}

std::unique_ptr<InteractionState> InteractionStateFactory::createInteractionState(Interaction *interaction,
                                                                                  std::string interactionStateType,
                                                                                  Cell *cell1,
                                                                                  Cell *cell2) {
    std::string identifier = interaction->getIdentifier();
    const auto&[type, next_states] = state_parameters_.at(identifier).at(interactionStateType);
    const auto clone = abm::util::overloaded{
            [&](InteractionType type) -> std::unique_ptr<InteractionType> {
                return std::make_unique<InteractionType>(type);
            },
            [&](Contacting type) -> std::unique_ptr<InteractionType> { return std::make_unique<Contacting>(type); },
            [&](RigidContacting type) -> std::unique_ptr<InteractionType> {
                return std::make_unique<RigidContacting>(type);
            },
            [&](Ingestion type) -> std::unique_ptr<InteractionType> { return std::make_unique<Ingestion>(type); },
    };

    std::unique_ptr<InteractionState>
            intState = std::make_unique<InteractionState>(interactionStateType, interaction, std::visit(clone, type),
                                                          next_states.empty());

    if (intState->getInteractionType() == "Ingestion") {
        intState = std::make_unique<InteractionState>(interactionStateType, interaction, std::visit(clone, type),
                                                      false);
        if (interaction->getFirstCell()->getTypeName() == "Macrophage") {
            interaction->getFirstCell()->addIngestions(interaction->getSecondCell()->getId());
        } else if (interaction->getSecondCell()->getTypeName() == "Macrophage") {
            interaction->getSecondCell()->addIngestions(interaction->getFirstCell()->getId());
        }
    }

    intState->addNextStateWithRate(next_states);

    return intState;
}

void InteractionStateFactory::addNextStates(Interaction *interaction,
                                            InteractionState *intState,
                                            const XMLNode &nextStates, Cell *cell1,
                                            Cell *cell2) {
    int noOfNextStates = nextStates.nChildNode();
    for (int i = 0; i < noOfNextStates; i++) {
        XMLNode curNextState = nextStates.getChildNode(i);
        std::string nameOfNextState = curNextState.getName();
        XMLNode rateOfNextStateNode;
        std::string curNextStateNextNode = curNextState.getChildNode(0).getName();
        if (curNextStateNextNode == "If") {
            XMLNode condition =
                    curNextState.getChildNode(0).getChildNode("Condition");
            for (int j = 0; j < condition.nChildNode(); j++) {
                if (cell1->getTypeName() == condition.getChildNode(j).getName()) {
                    InputConfiguration ic = InputConfiguration();
                    std::string nameOfConditon = condition.getChildNode(0).getChildNode(j).getName();
                    XMLNode conditionCell = condition.getChildNode(0).getChildNode(j);
                    int val = InputConfiguration::getIntDataFieldValueByName(&conditionCell, nameOfConditon);

                    if (cell1->getFeatureValueByName(nameOfConditon) == val) {
                        rateOfNextStateNode = curNextState.getChildNode(0).getChildNode("Then").getChildNode("Rate");
                        if (intState->getStateName() == "Phagocytose") {
                        }

                    } else {
                        rateOfNextStateNode =
                                curNextState.getChildNode(0).getChildNode("Else").getChildNode(
                                        "Rate");
                        if (intState->getStateName() == "Phagocytose") {
                        }
                    }
                } else if (cell2->getTypeName() == condition.getChildNode(j).getName()) {
                    std::string nameOfConditon = condition.getChildNode(0).getChildNode(j).getName();
                    XMLNode conditionCell = condition.getChildNode(0).getChildNode(j);
                    int val = InputConfiguration::getIntDataFieldValueByName(&conditionCell, nameOfConditon);
                    if (cell2->getFeatureValueByName(nameOfConditon) == val) {
                        rateOfNextStateNode = curNextState.getChildNode(0).getChildNode("Then").getChildNode("Rate");
                    } else {
                        rateOfNextStateNode =
                                curNextState.getChildNode(0).getChildNode("Else").getChildNode("Rate");
                    }
                }
            }

        } else {
            rateOfNextStateNode = curNextState.getChildNode("Rate");
        }
        std::string rateKey = InputConfiguration::getStringDataFieldByName(&rateOfNextStateNode, "ratekey");
        auto rate = RateFactory::getRate(rateKey);
        intState->addNextStateWithRate(nameOfNextState, rate);
    }
}
