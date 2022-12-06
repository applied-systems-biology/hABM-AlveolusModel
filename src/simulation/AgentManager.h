//  Copyright by Christoph Saffer, Paul Rudolph, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#ifndef AGENTMANAGER_H
#define AGENTMANAGER_H

#include <map>
#include <memory>
#include <set>
#include <utility>
#include <vector>

#include "io/XMLFile.h"
#include "simulation/morphology/SphereRepresentation.h"

class Site;
class Analyser;
class Agent;
class Cell;

class AgentManager {
public:
    /// Class for managing all the agents in data container. It provides functionality for input and replacement of agents during a simulation.
    AgentManager(double time_delta, Site *site);

    using iterator = std::vector<std::shared_ptr<Agent>>::iterator;
    using const_iterator = std::vector<std::shared_ptr<Agent>>::const_iterator;
    iterator begin() { return allAgents.begin(); }
    iterator end() { return allAgents.end(); }
    [[nodiscard]] const_iterator begin() const { return allAgents.begin(); }
    [[nodiscard]] const_iterator end() const { return allAgents.end(); }
    std::shared_ptr<Agent> &emplace_back(std::shared_ptr<Agent> &&value) {return allAgents.emplace_back(std::forward<std::shared_ptr<Agent>>(value));}

    /*!
     * Inputs an agent
     * @param current_time Double that contains current time
     * @param max_time Double that contains max time
     * @param random_generator Randomizer object that contains grandom generator
     */
    void inputOfAgents(double current_time, double max_time, Randomizer *random_generator);

    /*!
     * Inserts an agent at the systems boundary (e.g. alveolar entrance ring or pores of kohn)
     * @param site Site object (e.g. AlveoleSite)
     * @param agentType String that contains agent type
     */
    void insertAgentAtBoundary(Site *site, std::string agentType, double current_time);

    /*!
     * Replaces an agent with a new agent at a given position and move
     * @param site Site object (e.g. AlveoleSite)
     * @param agentToReplace Agent object that contains agent that is replaced
     * @param c Coordinate3D object that contains coordinates of new agent
     * @param prevMove Coordinate3D object that contains corrdinates of previous move
     */
    void replaceAgent(Site *site, Agent *agentToReplace, std::unique_ptr<Coordinate3D> c, Coordinate3D * prevMove, double current_time);

    /*!
     * Replaces an agent with another agent
     * @param site Site object (e.g. AlveoleSite)
     * @param agentToReplace Agent object that contains agent that is replaced
     * @param newAgent Agent object that replaces the previous agent
     */
    void replaceAgent(Site *site, Agent *agentToReplace, std::shared_ptr<Agent> newAgent, double current_time);

    /// Removes agent from the system
    void removeAgent(Site *site, Agent *agent, double current_time);

    /// Removes agents from the system that were previously set to deleted
    void cleanUpAgents();

    /// Removes a sphere object
    void removeSphereRepresentation(SphereRepresentation *sphereRep);
    void incrementIdHandling();
    void addConidiaToList(Agent *conidia) { posAfumiList.emplace_back(conidia); };
    void removeConidiaFromList(int id, double current_time);

    Agent *createAgent(Site *, std::string, Coordinate3D, double current_time);
    Agent *createAgent(Site *, std::string, Coordinate3D, Coordinate3D *, double current_time);
    int generateNewID() { return idHandling++;}

    void setLambdaInput(double speed, double persistenceTime);
    void setInitConQuantity() { initConQuantity = posAfumiList.size(); }
    void setLastConidiaChange(double lcc) { lastConidiaChange = lcc; }
    int getAgentQuantity(std::string agenttype);
    int getNextSphereRepresentationId(SphereRepresentation *sphereRep);
    [[nodiscard]] double getLastConidiaChange() const { return lastConidiaChange; };
    [[nodiscard]] int getIdHandling() const;
    [[nodiscard]] int getInitConQuantity() const { return initConQuantity; }
    double getOccupancyDensityOfSpace();
    Cell *getCellBySphereRepId(int sphereRepId);
    const std::vector<std::shared_ptr<Agent>> &getAllAgents();
    SphereRepresentation *getSphereRepBySphereRepId(int sphereRepId);
    std::set<SphereRepresentation *> *getAllSphereRepresentations() { return &allSphereRepresentations; };
    std::vector<Agent *> getAllConidia() { return posAfumiList; };
    std::vector<std::string> getAgentXMLTocTags(XMLFile *xml_file, bool csv) const;

private:
    std::vector<std::shared_ptr<Agent>> allAgents;
    std::map<int, Cell *> sphereIdToCell;
    std::map<int, SphereRepresentation *> sphereIdToSphereRep;
    std::set<SphereRepresentation *> allSphereRepresentations;
    int idHandling;
    int idHandlingSphereRepresentation;
    double lastConidiaChange{};
    double lastInputEventTime{};
    double nextInputEventTime{};
    unsigned int lastQuantity{};
    std::vector<Agent *> posAfumiList;
    int initConQuantity;
    double lambdaInput{};
    std::vector<double> conidiaRemoveTimes;
    std::vector<double> conidiaRemoveIDs;
    Site *site{};
    double time_delta_;
};

#endif    /* AGENTMANAGER_H */

