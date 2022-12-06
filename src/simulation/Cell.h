//  Copyright by Christoph Saffer, Paul Rudolph, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#ifndef CELL_H
#define    CELL_H

#include <memory>

#include "simulation/Agent.h"
#include "utils/io_util.h"
#include "simulation/Site.h"

class AgentProperties;
class Analyser;

class Cell : public Agent {
public:
    /// Class for an abstract cell that is further specified by inherited classes
    Cell(std::unique_ptr<Coordinate3D>, int, Site *, double time_delta, double current_time);

    /*!
     * Central function: Performs all actions for one timestep for one cell
     * @param timestep Double for courrent timestep
     * @param current_time Double for current time
     */
    void doAllActionsForTimestep(double timestep, double current_time) final;

    void setState(std::shared_ptr<CellState> cstate) final;
    void includeAgentXMLTagAoc(XMLFile *xmlTags, double current_time) final;
    void includeAgentXMLTagToc(XMLFile *xmlTags) final;
    void setPassive() override { passive = true; }
    void move(double timestep, double current_time) override;
    void setFeatureValueByName(std::string featureName, int value) override;
    void setVariableOnEvent(std::string variable, double value);
    void applyMethodByName(std::string mehtodName) override;
    void changeState(std::string stateName) override;
    void phagocytosisEvent(double time_delta, double current_time);
    void recieveInteractionEvent(InteractionEvent *ievent);
    void setInitialState(double time_delta, double current_time);
    void setExistingState(std::string stateName, double time_delta, double current_time);
    void printCellStates();
    void addIngestions(int id);

    int getIngestions() final { return ingestionCounter.size(); }
    double getFeatureValueByName(std::string featureName);
    Coordinate3D get_gradient() final { return cumulative_persistence_gradient; }
    Coordinate3D getEffectiveConnection(Cell *cell);
    std::string getTypeName() override;
    std::string getAgentCSVTagAoc(double current_time) final;
    std::string getAgentCSVTagToc() final;
    std::string generatePovObject() final;
    CellState *getCurrentCellState() final;
    std::shared_ptr<CellState> getCellStateByName(std::string nameOfState) final;
    size_t getIngestionPos(int id);
    Morphology *getSurface();
    Interactions *getInteractions();

    virtual void handleControlledAgents(double timestep);
    virtual void setCellInactive(double current_time);
    virtual void setPhagocytosed();
    virtual void passiveMove(double timestep, double current_time);
    virtual void doMorphologicalChanges(double timestep, double current_time);
    virtual void interactWithMolecules(double timestep);
    virtual void initialStateSetup(XMLNode cellNode);
    virtual bool cellActive();
    virtual void setup(double time_delta,
                       double current_time,
                       abm::util::SimulationParameters::AgentParameters *parameters);

protected:
    virtual void handleInteractionEvent(InteractionEvent *ievent);
    void setSurface();
    void setMovement(double timestep_size);
    void setPassiveMovement(double timestep_size);
    double phagocytosisEventDelay;
    std::shared_ptr<Morphology> surface;
    std::shared_ptr<Interactions> interactions;
    std::map<std::string, std::shared_ptr<CellState>> cellStates;
    std::shared_ptr<CellState> cellState;
    std::vector<int> ingestionCounter;

    //Molecule Parameters for AM
    std::map<std::string, double> secretion_rate;
    std::map<std::string, double> uptake_rate;
    std::map<std::string, double> k_blr;
    std::map<std::string, bool> chemotaxis;
    std::map<std::string, double> inside_conc_decay;
    std::map<std::string, double> k_r{};
    std::map<std::string, double> k_i{};
    double receptors{};
    double free_receptors{};
    double LR_complexes{};
    double R_internalized{};
    Coordinate3D cumulative_persistence_gradient;

};

#endif    /* CELL_H */

