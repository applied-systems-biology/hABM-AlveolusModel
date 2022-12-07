//  Copyright by Christoph Saffer, Paul Rudolph, Sandra Timme, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#include <numeric>
#include <algorithm>

#include "simulation/Site.h"
#include "simulation/Algorithms.h"
#include "simulation/boundary-condition/AbsorbingBoundaries.h"
#include "simulation/ParticleManager.h"
#include "simulation/AgentManager.h"
#include "utils/macros.h"
#include "Interaction.h"

Site::Site(double time_delta, unsigned int spatial_dimensions, Randomizer *random_generator,
           std::shared_ptr<InSituMeasurements> measurements) : dimensions(spatial_dimensions), random_generator_(
        random_generator), measurements_(std::move(measurements)), particle_manager_(
        std::make_unique<ParticleManager>(this)), agent_manager_(std::make_unique<AgentManager>(time_delta, this)) {

    measurements_->setSite(this);
}

void Site::setBoundaryCondition() {
    boundary_condition_ = std::make_unique<AbsorbingBoundaries>(this);
}

void Site::doAgentDynamics(Randomizer *random_generator, const SimulationTime &time) {

    const auto current_time = time.getCurrentTime();
    const auto dt = time.getCurrentDeltaT();
    const auto max_time = time.getMaxTime();

    const auto &all_agents = agent_manager_->getAllAgents();
    auto &all_particles = particle_manager_->getAllParticles();
    if (!all_agents.empty() || !all_particles.empty()) {
        // Loop over all agents (random order)
        const auto current_order = Algorithms::generateRandomPermutation(random_generator, all_agents.size());
        for (auto agent_idx = current_order.begin(); agent_idx < current_order.end(); ++agent_idx) {
            auto curr_agent = all_agents[*agent_idx];
            if (nullptr != curr_agent) {
                // Do all actions for one timestep for each agent (-> Cell.cpp)
                curr_agent->doAllActionsForTimestep(dt, current_time);
                // Remove spherical representations if the current agent got deleted
                if (curr_agent->isDeleted()) {
                    for (const auto &sphere: curr_agent->getAgentProperties()->getMorphology()->getAllSpheresOfThis()) {
                        neighbourhood_locator_->removeSphereRepresentation(sphere);
                        agent_manager_->removeSphereRepresentation(sphere);
                    }
                    curr_agent = nullptr;
                }
            }
        }

        // Loop over all particles if steady state is not reached
        if (!particle_manager_->steadyStateReached(current_time)) {
            all_particles = particle_manager_->getAllParticles();
            for (const auto &cur_particle : all_particles) {
                // Do all actions for one timestep for each particle
                cur_particle->doAllActionsForTimestep(dt);
            }
            // Initialize Particles
            particle_manager_->inputOfParticles(dt);
            for (const auto &cur_particle : all_particles) {
                // Apply actual concentration change to particles
                cur_particle->applyConcentrationChange(dt);
            }
        }

        // Clean up agents and particles
        agent_manager_->cleanUpAgents();
        agent_manager_->inputOfAgents(current_time, max_time, random_generator);
        measurements_->observeMeasurements(time);
    }
}

void Site::updateTimeStepSize(SimulationTime &time) {
    if (particle_manager_->getDiffusionCoefficient() > 500) {
        if (particle_manager_->steadyStateReached(time.getCurrentTime()) &&
            !large_timestep_active) {// && abs((systime_min) - round(systime_min)) < 0.000001) {
            // increase the timestep as particle dynamics are not affected anymore
            time.updateDeltaT(0.1);
            large_timestep_active = true;
            DEBUG_STDOUT("Increasing timestep to " << time.getCurrentDeltaT() << " at time " << time.getCurrentTime());
        }
        // If a conidia is found, reduce timestep again to process particle dynamics
        if (!particle_manager_->steadyStateReached(time.getCurrentTime()) && large_timestep_active) {
            time.updateDeltaT(time.getLastDeltaT());
            large_timestep_active = false;
            DEBUG_STDOUT("Decreasing timestep to " << time.getCurrentDeltaT() << " at time " << time.getCurrentTime());
        }
    }
}

bool Site::checkForStopping(const SimulationTime &time) const {
    return stopSimulation;
}

void Site::setStoppingCondition(const std::vector<std::string> &stopping_criteria) {
    for (const auto &criterion: stopping_criteria) {
        if (criterion == "FirstPassageTime") {
            stopping_FPT = true;
        }
    }
}

void Site::terminateSimulationForInteraction(Interaction &interaction) {

    // Ends a simulation if all fungi were touched at least once (FTP)
    // Cumulated first passage time (FTP) is clearance time (CT)
    if (stopping_FPT && interaction.getInteractionName() == "PhagocyteFungusInteraction") {
        int cellid{};
        if (interaction.getFirstCell()->getTypeName() == "AspergillusFumigatus") {
            cellid = cell_ids_FPT.emplace_back(interaction.getFirstCell()->getId());
        } else {
            cellid = cell_ids_FPT.emplace_back(interaction.getSecondCell()->getId());
        }
        cell_ids_FPT.emplace_back(cellid);
        DEBUG_STDOUT("Aspergillus Fumigatus with id=" << cellid << " was touched");
        std::sort(cell_ids_FPT.begin(), cell_ids_FPT.end());
        cell_ids_FPT.erase(std::unique(cell_ids_FPT.begin(), cell_ids_FPT.end()), cell_ids_FPT.end());

        int conidia_remaining = agent_manager_->getInitConQuantity() - cell_ids_FPT.size();
        if (conidia_remaining == 0 && stopping_FPT && agent_manager_->getInitConQuantity() > 0) {
            stopSimulation = true;
            DEBUG_STDOUT("All fungus were touched .. end simulation");
        }
    }
}
