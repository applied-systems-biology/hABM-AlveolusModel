//  Copyright by Christoph Saffer, Paul Rudolph, Sandra Timme, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#include <cmath>

#include "simulation/AgentManager.h"
#include "simulation/CellFactory.h"
#include "simulation/Interactions.h"
#include "simulation/cells/interaction/PhagocyteFungusInteraction.h"
#include "analyser/InSituMeasurements.h"
#include "utils/macros.h"
#include "simulation/RateFactory.h"
#include "io/InputConfiguration.h"
#include "simulation/Site.h"
#include "simulation/InteractionState.h"


AgentManager::AgentManager(double time_delta, Site *site) {
    this->site = site;
    time_delta_ = time_delta;
    idHandling = 0;
    idHandlingSphereRepresentation = 0;
}

void AgentManager::inputOfAgents(double current_time, double max_time, Randomizer *random_generator) {
    double lambda = site->getInputRate();
    if (current_time == 0) {
        // Generate the initial nextInputEventTime
        lastInputEventTime = 0;
        nextInputEventTime = 1.0 / (lambda) * log(1.0 / random_generator->generateDouble());
        DEBUG_STDOUT("Next AM Input:" + std::to_string(nextInputEventTime) + " with lambda of :" + std::to_string(lambda));
    }
    if (current_time >= nextInputEventTime) {
        //Generate nextInputEventTime when time the previously drawn nextInputEventTime is exceeded
        insertAgentAtBoundary(site, "Macrophage", current_time);
        lastInputEventTime = nextInputEventTime;
        double diff = 1.0 / (lambda) * log(1.0 / random_generator->generateDouble());
        nextInputEventTime = lastInputEventTime + diff;
        DEBUG_STDOUT("Time:" + std::to_string(current_time) + " - Next AM Input:" + std::to_string(nextInputEventTime));
    }
}

void AgentManager::insertAgentAtBoundary(Site *site, std::string agentType, double current_time) {
    Coordinate3D initialPosition;
    Coordinate3D initialVector;

    Agent *agent = nullptr;
    int rejections = 0;

    // Add an agent to the system without having any collisions
    do {
        if (agent != nullptr) {
            removeAgent(site, agent, current_time);
        }
        initialPosition = site->getRandomBoundaryPoint();
        initialVector = site->getBoundaryInputVector();
        agent = createAgent(site, agentType, initialPosition, &initialVector, current_time);
        rejections++;
    } while ((agent->getAgentProperties()->getInteractions()->hasCollisions() && rejections < 10000));
    if (rejections > 9999) {
        DEBUG_STDOUT("Could not find a position for the agent at the boundary. Agent is not added to the system.");
        if (agent != nullptr) removeAgent(site, agent, current_time);
    }
}

Agent *AgentManager::createAgent(Site *site,  std::string agenttype, Coordinate3D c,  double current_time) {
    auto agent = emplace_back(CellFactory::createCell(agenttype,  std::make_unique<Coordinate3D>(c),
                                                      generateNewID(),  site,  time_delta_, current_time));
    if (agenttype.compare("AspergillusFumigatus") == 0) {
        posAfumiList.push_back(agent.get());
    }
    return agent.get();
}

Agent *AgentManager::createAgent(Site *site, std::string agenttype, Coordinate3D c, Coordinate3D *prevMove, double current_time) {
    auto agent = emplace_back(CellFactory::createCell(agenttype, std::make_unique<Coordinate3D>(c), generateNewID(),
                                                      site, time_delta_, current_time));
    agent->getMovement()->setPreviousMove(prevMove);
    return agent.get();
}

void AgentManager::replaceAgent(Site *site,Agent *agent, std::unique_ptr<Coordinate3D> newCoord, Coordinate3D *prevMove,
                                double current_time) {

    std::string agentType = agent->getTypeName();
    auto newAgent = CellFactory::createCell(agentType, std::move(newCoord), idHandling, site, time_delta_, current_time);
    if (newAgent != 0) {
        newAgent->getMovement()->setPreviousMove(prevMove);
        std::replace_if(allAgents.begin(), allAgents.end(), [agent](const auto &a) { return agent == a.get(); },
                        newAgent);
        idHandling++;
    }
    agent->setDeleted();

}

void
AgentManager::replaceAgent(Site *site, Agent *agentToReplace, std::shared_ptr<Agent> newAgent, double current_time) {

    // If the new agent is the nullptr, delete the agent
    if (!newAgent) {
        auto all_interactions = agentToReplace->getAgentProperties()->getInteractions()->getAllInteractions();
        for (size_t i = 0; i < all_interactions.size(); i++) {
            if (all_interactions.at(i)->getInteractionName() == "PhagocyteFungusInteraction") {
                Cell *cell1 = all_interactions.at(i)->getFirstCell();
                Cell *cell2 = all_interactions.at(i)->getSecondCell();
                cell1->getTypeName() == "AspergillusFumigatus" ? cell1->setDeleted() : cell2->setDeleted();
            }
        }
        agentToReplace->setDeleted();
    } else {
        std::replace_if(allAgents.begin(),
                        allAgents.end(),
                        [agent = agentToReplace](const auto &a) { return agent == a.get(); },
                        newAgent);
    }

}

void AgentManager::removeAgent(Site *site, Agent *agent, double current_time) {

    // Remove agent and all its corresponding spheres in the neighbourhoodlocator
    site->getNeighbourhoodLocator()->removeSphereRepresentation(getSphereRepBySphereRepId(agent->getId()));
    agent->setDeleted();
    if (agent->getTypeName().compare("AspergillusFumigatus") == 0) {
        removeConidiaFromList(agent->getId(), current_time);
    }
    allAgents.erase(std::remove_if(allAgents.begin(),
                                   allAgents.end(),
                                   [agent](const auto &a) { return agent == a.get(); }), allAgents.end());
}

int AgentManager::getAgentQuantity(std::string agenttype) {
    int count = 0;
    for (auto agent: allAgents) {
        if (agent->getCurrentCellState()->getStateName() != "Death" && agent != 0) {
            if (agent->getTypeName() == agenttype) {
                count++;
            }
        }
    }
    return count;
}

const std::vector<std::shared_ptr<Agent>> &AgentManager::getAllAgents() {
    return allAgents;

}

void AgentManager::cleanUpAgents() {
    for (auto it = allAgents.begin(); it != allAgents.end();) {
        if (*it == 0) {
            it = allAgents.erase(it);
        } else {
            if ((*it)->isDeleted()) {
                for (const auto &sphere: (*it)->getAgentProperties()->getMorphology()->getAllSpheresOfThis()) {
                    site->getNeighbourhoodLocator()->removeSphereRepresentation(sphere);
                    removeSphereRepresentation(sphere);
                }
                it = allAgents.erase(it);
            } else {
                it++;
            }
        }
    }
}

int AgentManager::getNextSphereRepresentationId(SphereRepresentation *sphereRep) {
    sphereIdToCell[idHandlingSphereRepresentation] =
            sphereRep->getMorphologyElementThisBelongsTo()->getMorphologyThisBelongsTo()->getCellThisBelongsTo();
    sphereIdToSphereRep[idHandlingSphereRepresentation] = sphereRep;
    allSphereRepresentations.insert(sphereRep);
    return idHandlingSphereRepresentation++;
}

Cell *AgentManager::getCellBySphereRepId(int sphereRepId) {
    if (sphereIdToCell.find(sphereRepId) != sphereIdToCell.end()) {
        return sphereIdToCell[sphereRepId];
    }
    return nullptr;
}

SphereRepresentation *AgentManager::getSphereRepBySphereRepId(int sphereRepId) {
    if (sphereIdToSphereRep.find(sphereRepId) != sphereIdToSphereRep.end()) {
        return sphereIdToSphereRep[sphereRepId];
    }
    return nullptr;
}

void AgentManager::removeSphereRepresentation(SphereRepresentation *sphereRep) {
    allSphereRepresentations.erase(sphereRep);
    sphereIdToSphereRep.erase(sphereRep->getId());
    sphereIdToCell.erase(sphereRep->getId());
}

void AgentManager::setLambdaInput(double speed, double persistenceTime) {
    double lambda;
    if (lambdaInput == 0) {
        std::ostringstream ssSpeed, ssTp;
        ssSpeed << speed;
        ssTp << persistenceTime;
        InputConfiguration ic = InputConfiguration("input/parameter/ED-lambda.xml");
        XMLNode lambdas = ic.getRootNode().getChildNode("lambda-input-values");
        const char *lambdaChar = lambdas
                        .getChildNodeWithAttribute("speed", "value",
                                                   ssSpeed.str().c_str())
                        .getChildNodeWithAttribute("tp", "value", ssTp.str().c_str())
                        .getChildNode("lambda")
                        .getAttribute("value");
        lambda = atof(lambdaChar);
        lambdaInput = lambda;
    }
}

int AgentManager::getIdHandling() const {
    return idHandling;
}

void AgentManager::incrementIdHandling() {
    idHandling++;
}

double AgentManager::getOccupancyDensityOfSpace() {
    double cellVolume = 0;
    Coordinate3D lowerLim = site->getLowerLimits();
    Coordinate3D upperLim = site->getUpperLimits();

    double x = abs(lowerLim.x) + upperLim.x;
    double y = abs(lowerLim.y) + upperLim.y;
    double z = abs(lowerLim.z) + upperLim.z;
    double siteVolume = 0;
    if (site->getNumberOfSpatialDimensions() == 2) {
        siteVolume = x * y;
    } else if (site->getNumberOfSpatialDimensions() == 3) {
        siteVolume = x * y * z;
    }

    for (size_t i = 0; i < allAgents.size(); i++) {
        cellVolume = cellVolume + allAgents[i]->getAgentProperties()->getMorphology()->getVolume();
    }
    double occupancy = (cellVolume * 100) / siteVolume;
    return occupancy;
}

void AgentManager::removeConidiaFromList(int id, double current_time) {
    for (size_t i = 0; i < posAfumiList.size(); i++) {
        if (posAfumiList.at(i)->getId() == id) {
            posAfumiList.erase(posAfumiList.begin() + i);
            conidiaRemoveTimes.push_back(current_time);
            conidiaRemoveIDs.push_back(id);
            // When a conidia is found current time must be stored to decide about steady state
            setLastConidiaChange(current_time);
            // Particles must be updated
            DEBUG_STDOUT("Conidia removed. Possible Chemotaxis cleanup initiated at " << current_time);
            site->getParticleManager()->setCleanChemotaxis(true);
            break;
        }
    }
}

std::vector<std::string> AgentManager::getAgentXMLTocTags(XMLFile *xml_node, bool csv) const {
    std::vector<std::string> csv_tags;
    if (csv) csv_tags.reserve(allAgents.size());
    for (const auto &agent: allAgents) {
        agent->includeAgentXMLTagToc(xml_node);
        if (csv) csv_tags.emplace_back(agent->getAgentCSVTagToc());
    }
    return csv_tags;
}