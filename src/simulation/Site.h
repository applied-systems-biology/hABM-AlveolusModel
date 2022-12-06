//  Copyright by Christoph Saffer, Paul Rudolph, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#ifndef _SITE_H
#define _SITE_H

#include <memory>
#include <utility>
#include <vector>
#include <algorithm>

#include "utils/time_util.h"
#include "basic/Coordinate3D.h"
#include "simulation/Cell.h"
#include "simulation/Algorithms.h"
#include "io/InputConfiguration.h"
#include "analyser/InSituMeasurements.h"
#include "simulation/morphology/SphereRepresentation.h"
#include "simulation/neighbourhood/BoundaryCondition.h"
#include "simulation/NeighbourhoodLocator.h"
#include "simulation/ParticleManager.h"
#include "simulation/AgentManager.h"
#include "io/output_handler.h"

class Agent; //forward declaration
class Particle; //forward declaration
class Analyser;
class Randomizer;

class Site {

public:
    /// Class for handling environment interactions and wrapping functionality of all aspects that are happening during the simulation inside of the environment (i.e. site)
    Site(double time_delta,
         unsigned int spatial_dimensions,
         Randomizer *random_generator,
         std::shared_ptr<InSituMeasurements> measurements);

    virtual ~Site() = default;

    /*!
     * Conducts all actions that happen in one timestep in the site
     * @param random_generator Randomizer object that contains randomizer of current run
     * @param time SimulationTime object, i.e. contains current time and timestep
     */
    void doAgentDynamics(Randomizer *random_generator, const SimulationTime &time);

    /*!
     * Updates the timestep size (i.e. steady state of diffusion was reached)
     * @param time SimulationTime object, i.e. contains current time and timestep
     */
    void updateTimeStepSize(SimulationTime &time);

    /*!
     * Sets stopping conditions of simulator-config.json
     * @param stopping_criteria Vector of Strings of stopping criteria (i.e. "stopping_criteria": ["FirstPassageTime"])
     */
    void setStoppingCondition(const std::vector<std::string> &stopping_criteria);

    /*!
     * Terminates Simulation for certain interactions (i.e. all phagocytes were touched)
     * @param interaction Interaction object that contains an interaction
     */
    void terminateSimulationForInteraction(Interaction &interaction);

    [[nodiscard]] bool checkForStopping(const SimulationTime &time) const;
    Randomizer *getRandomGenerator() { return random_generator_; }
    NeighbourhoodLocator *getNeighbourhoodLocator() { return neighbourhood_locator_.get(); }
    ParticleManager *getParticleManager() const { return particle_manager_.get(); }
    AgentManager *getAgentManager() const { return agent_manager_.get(); }
    InSituMeasurements *getMeasurments() const { return measurements_.get(); }
    Coordinate3D getBoundaryInputVector() { return boundary_input_vector_; }
    bool getLargeTimestepActive() { return large_timestep_active; }
    [[nodiscard]] int getState() const { return state_; }
    [[nodiscard]] unsigned int getNumberOfSpatialDimensions() const { return dimensions; }
    [[nodiscard]] int getBoundaryInput() const { return boundaryInput; }
    [[nodiscard]] int getBoundaryParticleInput() const { return boundaryParticleInput; }
    [[nodiscard]] double getLatestAlpha2dTurningAngle() const { return alpha2dTurningAngle; }
    [[nodiscard]] double getInputRate() const { return inputRate; }
    [[nodiscard]] std::string getIdentifier() const { return identifier_; }

    friend void OutputHandler::outputCurrentConfiguration(const Site &site,
                                                          const SimulationTime &time,
                                                          int run,
                                                          int seed,
                                                          bool simulation_end,
                                                          std::unordered_map<std::string, std::string> cmd_input_args) const;
    friend void InSituMeasurements::observeMeasurements(const SimulationTime &time);

    virtual void handleBoundaryCross(Agent *, Coordinate3D *, double current_time) = 0;
    virtual void includeSiteXMLTagToc(XMLFile *xml_file) const = 0;
    virtual bool containsPosition(Coordinate3D position) = 0;
    [[nodiscard]] virtual std::string getType() const = 0;
    virtual Coordinate3D getRandomPosition() = 0;
    virtual Coordinate3D getCenterPosition() = 0;
    virtual Coordinate3D getRandomBoundaryPoint() = 0;
    virtual Coordinate3D getLowerLimits() = 0;
    virtual Coordinate3D getUpperLimits() = 0;
    virtual Coordinate3D generateRandomDirectionVector(Coordinate3D position, double length) = 0;
    virtual Coordinate3D generatePersistentDirectionVector(Coordinate3D position,
                                                           double length,
                                                           Coordinate3D prevVector,
                                                           double previousAlpha) = 0;
    virtual Coordinate3D generateBackShiftOnContacting(SphereRepresentation *activeSphere,
                                                       SphereRepresentation *passiveSphere,
                                                       double mustOverhead) = 0;
    virtual void dynamicChange(double timestep, double current_time) {};
    virtual bool overAECT1(SphericCoordinate3D posConida) { return {}; }
    virtual bool onAECTObstacleCell(Coordinate3D position) { return {}; };
    virtual double getFeatureValueByName(std::string name) { return {}; }
    [[nodiscard]] virtual double getRadius() const { return 0.0; }
    virtual double getThicknessOfBorder() { return 0.0; }
    virtual double getLowerThetaBound() { return 0.0; }
    virtual double getDistanceFromBoundary(Coordinate3D position) { return 0.0; }
    virtual Coordinate3D getRandomMinDistanceToBoundaryPosition(double minDistanceToBoundary) { return {}; }
    virtual Coordinate3D generateBiasedRandomDirectionVector(Agent *agent,
                                                             Coordinate3D position,
                                                             double length) { return {}; }
    virtual Coordinate3D generateDirectedVector(Coordinate3D position,
                                                SphericCoordinate3D posOfGoal,
                                                double length) { return {}; }
    virtual Coordinate3D generateDirectedVector(Coordinate3D position,
                                                double alpha, double length) { return {}; }
    virtual std::vector<SphericCoordinate3D> getAECT1() { return {}; };
    virtual std::vector<SphericCoordinate3D> getAECT2() { return {}; };
    virtual std::vector<SphericCoordinate3D> getPOK() { return {}; };
    virtual std::vector<Coordinate3D> getSystemBoundaries() { return {}; };

protected:
    void setBoundaryCondition();
    virtual void initializeAgents(const abm::util::SimulationParameters::AgentManagerParameters &parameters,
                          const std::string &input_dir,
                          double current_time,
                          double time_delta){};
    int state_{};
    int boundaryInput{};
    int boundaryParticleInput{};
    unsigned int dimensions{};
    bool passiveMovementOn{};
    double inputRate{};
    double alpha2dTurningAngle{};
    std::string identifier_{};
    bool large_timestep_active = false;
    bool stopSimulation{};
    bool stopping_FPT = false;
    std::vector<int> cell_ids_FPT{};
    Coordinate3D boundary_input_vector_{};
    std::vector<std::pair<std::string, long>> stopping_cell_states;
    Randomizer *random_generator_;
    std::shared_ptr<InSituMeasurements> measurements_;
    std::unique_ptr<BoundaryCondition> boundary_condition_;
    std::unique_ptr<NeighbourhoodLocator> neighbourhood_locator_;
    std::unique_ptr<ParticleManager> particle_manager_;
    std::unique_ptr<AgentManager> agent_manager_;
};

#endif    /* _SITE_H */

