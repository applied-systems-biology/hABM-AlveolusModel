//  Copyright by Christoph Saffer, Paul Rudolph, Sandra Timme, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#ifndef _AGENT_H
#define    _AGENT_H

#include <iostream>

#include "basic/Coordinate3D.h"
#include "simulation/movement/RandomWalk.h"
#include "simulation/Movement.h"
#include "simulation/Morphology.h"
#include "simulation/cells/CellState.h"
#include "io/XMLFile.h"
#include "simulation/AgentProperties.h"

class Site; //forward declaration
class Interactions; //forward declaration
class HyphalBranch; //forward declaration

class Agent {
public:
  // Abstract class for agents in the hABM. This class provides the cell class with it main functionality.
    Agent();
    Agent(std::unique_ptr<Coordinate3D>, int, Site *);
    virtual ~Agent() = default;

    void setPosition(Coordinate3D newPos);
    void setBeenMovedThisTimestep(bool newHasBeenMoved);
    void setTimestepLastTreatment(double current_time);
    void setDeleted();
    void setPreviousPosition(Coordinate3D *pos) { *previousPosition = *pos; }
    void resetAgent(Coordinate3D, double current_time);

    [[nodiscard]] int getId() const;
    [[nodiscard]] bool isDeleted() const { return is_deleted_; }
    [[nodiscard]] bool agentTreatedInCurrentTimestep(double current_time) const;
    bool hasBeenMovedThisTimestep();
    bool coordinateIsInsideAgent(Coordinate3D *, Agent *);
    bool hasCollisionsInsideAgent(Agent *, const std::vector<Agent *> &);
    bool shiftPosition(Coordinate3D *shifter,
                       double current_time,
                       SphereRepresentation *sphereRep = 0,
                       std::string originCall = "not known");
    double getLifetime(double);
    double getInitialTime();
    Site *getSite();
    AgentProperties *getAgentProperties();
    Movement *getMovement();
    Movement *getPassiveMovement();
    Coordinate3D getPosition();
    Coordinate3D *getCurrentShift();
    Coordinate3D getPreviousPosition() { return *previousPosition; };
    Coordinate3D getInitialPosition() { return *initialPosition; };
    Coordinate3D getCurrentPosition() { return *position; };
    Coordinate3D getCoordinateWithinAgent(Agent *);
    std::map<std::string, double> molecule_uptake;

    virtual void doAllActionsForTimestep(double timestep, double current_time) = 0;
    virtual void move(double timestep, double current_time) = 0;
    virtual std::string getTypeName() = 0;
    virtual CellState *getCurrentCellState() = 0;
    virtual std::string generatePovObject() = 0;
    virtual void setPassive() = 0;
    virtual void includeAgentXMLTagAoc(XMLFile *xmlTags, double current_time) = 0;
    virtual void includeAgentXMLTagToc(XMLFile *xmlTags) = 0;
    virtual std::string getAgentCSVTagAoc(double current_time) = 0;
    virtual std::string getAgentCSVTagToc() = 0;
    virtual void setVariableOnEvent(std::string variable, double value) = 0;
    virtual void setFeatureValueByName(std::string featureName, int value) = 0;
    virtual void applyMethodByName(std::string mehtodName) = 0;
    virtual void changeState(std::string stateName) = 0;
    virtual int getIngestions() = 0;
    virtual double getFeatureValueByName(std::string featureName) = 0;
    virtual std::shared_ptr<CellState> getCellStateByName(std::string nameOfState) = 0;
    virtual void setState(std::shared_ptr<CellState> state) = 0;
    virtual Morphology *getSurface() = 0;
    virtual Coordinate3D get_gradient() = 0;

protected:
    void setInitialPosition(Coordinate3D initPos);
    unsigned int id{};
    bool positionShiftAllowed;
    bool hasBeenMoved;
    bool passive;
    bool PoKset;
    bool is_deleted_;
    double initialTime;
    double timestepLastTreatment;

    Site *site;
    std::unique_ptr<AgentProperties> agentProps;
    std::unique_ptr<Coordinate3D> currShift;
    std::unique_ptr<Coordinate3D> initialPosition;
    std::unique_ptr<Coordinate3D> previousPosition;
    std::shared_ptr<Coordinate3D> position;
    std::shared_ptr<Movement> movement;
    std::shared_ptr<Movement> passiveMovement;
};

#endif    /* _AGENT_H */

