//  Copyright by Christoph Saffer, Paul Rudolph, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#include "simulation/Cell.h"

#include "simulation/Agent.h"
#include "simulation/AgentProperties.h"
#include "analyser/Analyser.h"
#include "simulation/Interaction.h"
#include "simulation/InteractionState.h"
#include "simulation/Interactions.h"
#include "simulation/CellStateFactory.h"
#include "simulation/morphology/SphericalMorphology.h"
#include "simulation/movement/BiasedPersistentRandomWalk.h"
#include "analyser/InSituMeasurements.h"
#include "utils/macros.h"


Cell::Cell(std::unique_ptr<Coordinate3D> c, int id, Site *site, double time_delta, double current_time)
        : Agent(std::move(c), id, site) {

    cellState = 0;
    interactions = std::make_shared<Interactions>(this, site->getNeighbourhoodLocator());
    agentProps->setInteractions(interactions);
    phagocytosisEventDelay = 0;
}

void Cell::doAllActionsForTimestep(double timestep, double current_time) {
    enum class Task {
        MOVEMENT,
        INTERACTIONS_AND_STATES,
        MORPHOLOGY_CHANGE,
        MOLECULE_INTERACTION
    };
    enum class Change {
        STATES,
        INTERACTIONS
    };

    // All actions for one cell in one timestep
    if (agentTreatedInCurrentTimestep(current_time)) {
        if (!is_deleted_) move(timestep, current_time);
    } else {
        std::vector<unsigned int>::iterator it;
        std::vector<unsigned int> currPerm = Algorithms::generateRandomPermutation(site->getRandomGenerator(), 4);
        unsigned currentTask;

        // Randomly ordered execution of Movement, Interaction_and_States, etc. per cell per timestep
        for (it = currPerm.begin(); it < currPerm.end(); it++) {
            currentTask = *it;
            if (is_deleted_) break;
            switch (static_cast<Task>(currentTask)) {
                case Task::MOVEMENT:
                    if (!is_deleted_) {
                        // do movement in current timestep of current cell
                        move(timestep, current_time);
                        if (is_deleted_) it = currPerm.end();
                    }
                    break;
                case Task::INTERACTIONS_AND_STATES: {
                    std::vector<unsigned int>::iterator statesIt;
                    std::vector<unsigned int> states = Algorithms::generateRandomPermutation(site->getRandomGenerator(),2);
                    for (statesIt = states.begin(); statesIt < states.end(); statesIt++) {
                        switch (static_cast<Change>(*statesIt)) {
                            case Change::STATES:
                                if (!agentTreatedInCurrentTimestep(current_time) && !is_deleted_) {
                                    // do state transiation in current timestep of current cell
                                    cellState->stateTransition(timestep, current_time);
                                }
                                break;
                            case Change::INTERACTIONS:
                                if (!agentTreatedInCurrentTimestep(current_time) && !is_deleted_) {
                                    // do all interactions in current timestep of current cell
                                    interactions->doWholeProcess(timestep, current_time, site->getMeasurments());
                                }
                                break;
                        }
                    }
                }
                    break;
                case Task::MORPHOLOGY_CHANGE:
                    // Morphology changes in current timestep of current cell (e.g. swelling of fungus)
                    doMorphologicalChanges(timestep, current_time);
                    break;
                case Task::MOLECULE_INTERACTION:
                    // Interaction with molecules in current timestep of current cell (e.g. AM)
                    interactWithMolecules(timestep);
                    break;
            }
        }
    }
}

void Cell::move(double timestep, double current_time) {
    Coordinate3D *move = movement->move(timestep, -1);
    shiftPosition(move, current_time, 0, "Movement");
    if (!is_deleted_ && move->getMagnitude() > 0) {
        interactions->avoidNewInteractions(timestep, current_time);
    }
}

void Cell::passiveMove(double timestep, double current_time) {
    Coordinate3D *pMove = passiveMovement->move(timestep, -1);
    shiftPosition(pMove, current_time, 0, "Movement");
    if (!is_deleted_ && pMove->getMagnitude() > 0) {
        interactions->avoidNewInteractions(timestep, current_time);
    }
}

void Cell::doMorphologicalChanges(double timestep, double current_time) {
}

std::string Cell::getTypeName() {
    return "Cell";
}

void Cell::setMovement(double timestep_size) {
}

void Cell::setPassiveMovement(double timestep_size) {
}


Interactions *Cell::getInteractions() {
    return interactions.get();
}

void Cell::setSurface() {
}

std::string Cell::generatePovObject() {
    return surface->generatePovObject();
}

void Cell::handleControlledAgents(double timestep) {
}

void Cell::phagocytosisEvent(double time_delta, double current_time) {
    phagocytosisEventDelay = current_time + time_delta * 50;
}

void Cell::includeAgentXMLTagAoc(XMLFile *xmlTags, double current_time) {
    std::ostringstream ssTime;
    XMLNode agentNode = xmlTags->addChildToRootNode("Agent");
    xmlTags->addDataFieldToNode(agentNode, *position);
    ssTime << current_time;
    xmlTags->addDataFieldToNode(agentNode, "time", "discrete", "double", ssTime.str());
    xmlTags->addDataFieldToNode(agentNode, "typename", "discrete", "string", getTypeName());
    XMLNode cellNode = xmlTags->addChildToNode(agentNode, "Cell");
    xmlTags->addDataFieldToNode(cellNode, "CellState", "discrete", "string", cellState->getStateName());
    surface->includeMorphologyXMLOutput(xmlTags, &cellNode);
    interactions->includeInteractionsXMLOutput(xmlTags, &cellNode);
}

void Cell::includeAgentXMLTagToc(XMLFile *xmlTags) {

    std::ostringstream ssid;
    XMLNode agentNode = xmlTags->addChildToRootNode("Agent");
    xmlTags->addDataFieldToNode(agentNode, *position);
    ssid << id;
    xmlTags->addDataFieldToNode(agentNode, "id", "discrete", "unsigned int", ssid.str());
    xmlTags->addDataFieldToNode(agentNode, "typename", "discrete", "string", getTypeName());
    XMLNode cellNode = xmlTags->addChildToNode(agentNode, "Cell");
    xmlTags->addDataFieldToNode(cellNode, "state", "discrete", "string", cellState->getStateName());
    surface->includeMorphologyXMLOutput(xmlTags, &cellNode);
    interactions->includeInteractionsXMLOutput(xmlTags, &cellNode);
}

std::string Cell::getAgentCSVTagAoc(double current_time) {
    std::ostringstream csvTag;
    csvTag << current_time << "\t" <<
           position->x << "\t" <<
           position->y << "\t" <<
           position->z << "\t" <<
           getTypeName();
    return csvTag.str();
}

std::string Cell::getAgentCSVTagToc() {

    std::ostringstream csvTag;
    csvTag << id << "\t" <<
           position->x << "\t" <<
           position->y << "\t" <<
           position->z << "\t" <<
           getTypeName();
    return csvTag.str();
}

void Cell::setState(std::shared_ptr<CellState> cstate) {
    cellState = cstate;
}

Coordinate3D Cell::getEffectiveConnection(Cell *cell) {
    return Coordinate3D(*cell->position) - Coordinate3D(*position);
}

void Cell::recieveInteractionEvent(InteractionEvent *ievent) {
    cellState->handleInteractionEvent(ievent);
    handleInteractionEvent(ievent);
}

void Cell::handleInteractionEvent(InteractionEvent *ievent) {
}

void Cell::initialStateSetup(XMLNode cellNode) {
}

void Cell::setInitialState(double time_delta, double current_time) {
}

void Cell::setExistingState(std::string stateName, double time_delta, double current_time) {
    if (cellStates.find(stateName) == cellStates.end()) {
        cellState = CellStateFactory::createCellState(this, stateName);
    } else {
        cellState = cellStates[stateName];
    }
    cellState->stateTransition(time_delta, current_time);
}

std::shared_ptr<CellState> Cell::getCellStateByName(std::string nameOfState) {
    if (cellStates.find(nameOfState) == cellStates.end()) {
        return nullptr;
    } else {
        return cellStates[nameOfState];
    }
}

CellState *Cell::getCurrentCellState() {
    return cellState.get();
}

double Cell::getFeatureValueByName(std::string featureName) {
    double value = 0;
    if (featureName.compare("migration-bias-probability") == 0) { // =p
        double LRdiff_front_rear = cumulative_persistence_gradient.getMagnitude();
        // According to the paper, should something like: value = LRdiff_front_rear * 1.2e-3
        value = LRdiff_front_rear * 0.6 / (0.424 * 500.0); // alternative proportionality factor calculated from Farrell et al. (1990)
        if (value > 1.0) {
            value = 1.0;
        }
    }

    if (featureName.compare("receptor-difference-front-rear") == 0) {
        value = cumulative_persistence_gradient.getMagnitude();
    }
    return value;
}

void Cell::printCellStates() {
    auto it = cellStates.begin();
    DEBUG_STDOUT("-------[CellStates]-------");
    while (it != cellStates.end()) {
        DEBUG_STDOUT(it->first);
        it++;
    }
    DEBUG_STDOUT("--------------------------");
}

size_t Cell::getIngestionPos(int id) {
    for (size_t i = 0; i < ingestionCounter.size(); i++) {
        if (ingestionCounter.at(i) == id) {
            return i;
            break;
        }
    }
    return 0;
}

void Cell::setCellInactive(double current_time) {
}

void Cell::setPhagocytosed() {
}

void Cell::interactWithMolecules(double timestep) {
}

bool Cell::cellActive() {
    return false;
}

void Cell::addIngestions(int id) {
    if (std::find(ingestionCounter.begin(), ingestionCounter.end(), id) != ingestionCounter.end()) {
        /* v already inside */
    } else {
        /* v not yet inside */
        ingestionCounter.push_back(id);
    }
}

Morphology *Cell::getSurface() {
    return surface.get();
}

void Cell::changeState(std::string stateName) {
}

void Cell::applyMethodByName(std::string mehtodName) {
}

void Cell::setVariableOnEvent(std::string variable, double value) {
    if (variable.compare("reset-cumulative-gradient") == 0) {
        cumulative_persistence_gradient *= 0.0;
    }
}

void Cell::setFeatureValueByName(std::string featureName, int value) {
}

void Cell::setup(double time_delta, double current_time, abm::util::SimulationParameters::AgentParameters *parameters) {

    if (parameters->movement_parameters.type.empty()) {
        movement = std::make_unique<Movement>(site->getNumberOfSpatialDimensions());
    } else if (parameters->movement_parameters.type == "RandomWalk") {
        if (parameters->movement_parameters.mean == 0) {
            auto speed = sqrt(
                    (2 * site->getNumberOfSpatialDimensions() * parameters->movement_parameters.diffusion_coefficient) /
                    time_delta);
            movement = std::make_unique<RandomWalk>(site, speed, 0, site->getNumberOfSpatialDimensions());
        } else {
            movement = std::make_unique<RandomWalk>(site,
                                                    parameters->movement_parameters.mean,
                                                    parameters->movement_parameters.stddev,
                                                    site->getNumberOfSpatialDimensions());
        }
    } else if (parameters->movement_parameters.type == "BiasedPersistentRandomWalk") {
        movement = std::make_unique<BiasedPersistentRandomWalk>(this,
                                                                parameters->movement_parameters.persistence_time,
                                                                parameters->movement_parameters.mean,
                                                                site->getNumberOfSpatialDimensions());
    }
    movement->setSite(site);
    movement->setCurrentPosition(position.get());
    agentProps->setMovement(movement);
    if (parameters->passive_movement_parameters.type.empty()) {
        passiveMovement = std::make_unique<Movement>(site->getNumberOfSpatialDimensions());
    } else if (parameters->passive_movement_parameters.type == "RandomWalk") {
        if (parameters->passive_movement_parameters.mean == 0) {
            auto speed =
                    sqrt((2 * site->getNumberOfSpatialDimensions() *
                          parameters->passive_movement_parameters.diffusion_coefficient) / time_delta);
            passiveMovement = std::make_unique<RandomWalk>(site, speed, 0);
        } else {
            passiveMovement =
                    std::make_unique<RandomWalk>(site, parameters->passive_movement_parameters.mean,
                                                 parameters->passive_movement_parameters.stddev);
        }
    } else if (parameters->passive_movement_parameters.type == "BiasedPersistentRandomWalk") {
        passiveMovement = std::make_unique<BiasedPersistentRandomWalk>(this,
                                                                       parameters->passive_movement_parameters.persistence_time,
                                                                       parameters->passive_movement_parameters.mean,
                                                                       site->getNumberOfSpatialDimensions());
    }

    passiveMovement->setSite(site);
    passiveMovement->setCurrentPosition(position.get());
    agentProps->setPassiveMovement(passiveMovement);
    timestepLastTreatment = -1;
    std::string color = parameters->morphology_parameters.color;
    surface = std::make_shared<Morphology>(color, this);
    if (parameters->morphology_parameters.type == "SphericalMorphology") {
        auto radius = site->getRandomGenerator()->generateNormalDistributedValue(
                parameters->morphology_parameters.radius,
                parameters->morphology_parameters.stddev);
        surface->appendAssociatedCellpart(std::make_unique<SphericalMorphology>(surface.get(),
                                                                                position,
                                                                                radius,
                                                                                "basic"));
    }
    agentProps->setMorphology(surface);

    for (const auto&[state, next_states]:CellStateFactory::getStateSetup(site->getIdentifier(), parameters->type)) {
        cellStates[state] = std::make_shared<CellState>(state, this, next_states);
    }
    if (cellStates.find("InitialCellState") == cellStates.end()) {
        cellState = CellStateFactory::createCellState(this, "InitialCellState");
    } else {
        cellState = cellStates["InitialCellState"];
    }
    cellState->stateTransition(time_delta, current_time);
    for (auto &[name, data]:parameters->molecule_interactions) {
        if (data.secretion != 0) {
            secretion_rate.emplace(std::pair<std::string, double>(name, data.secretion));
        }
        if (data.uptake != 0) {
            uptake_rate.emplace(std::pair<std::string, double>(name, data.uptake));
        }
        if (data.k_blr != 0) {
            k_blr.emplace(std::pair<std::string, double>(name, data.k_blr));
        }
        if (data.chemotaxis != 0) {
            chemotaxis.emplace(std::pair<std::string, bool>(name, data.chemotaxis));
        }
        if (data.k_i != 0) {
            k_i.emplace(std::pair<std::string, double>(name, data.k_i));
        }
        if (data.k_r != 0) {
            k_r.emplace(std::pair<std::string, double>(name, data.k_r));
        }
    }
}
