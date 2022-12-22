//  Copyright by Christoph Saffer, Paul Rudolph, Sandra Timme, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#include <chrono>
#include <omp.h>
#include <string>

#include "simulation/simulator.h"
#include "io/output_handler.h"
#include "simulation/CellFactory.h"
#include "simulation/InteractionFactory.h"
#include "simulation/RateFactory.h"
#include "utils/macros.h"
#include "utils/misc_util.h"
#include "utils/time_util.h"
#include "visualisation/Visualizer.h"
#include "simulation/Site.h"
#include "simulation/site/AlveoleSite.h"
#include "simulation/site/SphereSite.h"
#include "simulation/AgentManager.h"
#include "simulation/CellStateFactory.h"
#include "simulation/InteractionStateFactory.h"

int Simulator::consumers = 0;

Simulator::Simulator(std::string config_path, const std::unordered_map<std::string, std::string> &cmd_input_args)
        : config_path_(config_path) {

    // Read in configuration files and parameters
    const auto simulator_config = static_cast<boost::filesystem::path>(config_path).append("simulator-config.json");
    const auto input_config = static_cast<boost::filesystem::path>(config_path).append("input-config.json");
    parameters_ = abm::util::getSimulationParameters(simulator_config.string());
    parameters_.cmd_input_args = cmd_input_args;
    if (cmd_input_args.size() > 0) handleCmdInputArgs(cmd_input_args);
    auto input_parameters = abm::util::getInputParameters(input_config.string());
    if (parameters_.use_interactions) {
        InteractionFactory::initialize(parameters_.interaction_parameters);
    }
    // Prepare simulations
    RateFactory::initialize(input_parameters.rates);
    CellFactory::initialize(parameters_.site_parameters);
    CellStateFactory::initialize(parameters_.site_parameters);
    InteractionStateFactory::initialize(parameters_.interaction_parameters);
    ++Simulator::consumers;
}

Simulator::~Simulator() {
    --Simulator::consumers;
    if (parameters_.use_interactions) {
        InteractionFactory::close();
        InteractionStateFactory::close();
    }
    RateFactory::close();
    CellStateFactory::close();
    CellFactory::close();
}

void Simulator::executeRuns(int runs, int seed, const std::string &output_dir, const std::string &input_dir, int sim,
                            const std::string &parameter_string) const {
    std::ostringstream project_name;

    // Initializes the seed for the current parameter configuration
    SYSTEM_STDOUT("System Seed: " << seed);
    int current_sim_seed = seed + runs * sim;
    SYSTEM_STDOUT("Current Simulation Seed: " << current_sim_seed);

    // Initializes project name and output directory
    project_name << abm::util::getCurrentLocalTimeAsString() << "_" << parameter_string << current_sim_seed;
    const auto project_dir = static_cast<boost::filesystem::path>(output_dir).append(parameters_.topic).append(
            project_name.str()).string();

    // Initializes output handler, visualizer, analyzer and timer
    const auto output_handler = std::make_unique<const OutputHandler>(config_path_, project_dir);
    const auto visualizer = std::make_unique<const Visualizer>(config_path_, project_dir, runs);
    const auto analyser = std::make_unique<const Analyser>(config_path_, project_dir);
    using timer = std::chrono::steady_clock;
    const auto start = timer::now();

    // Start parallelized for-loop over all runs for one parameter configuration
#pragma omp parallel for schedule(dynamic)
    for (int current_run = 1; current_run <= runs; ++current_run) {
        SYSTEM_STDOUT("Thread " << omp_get_thread_num() << ": Start Run " << current_run << "/" << runs);

        // Setup environment for each run, e.g. each run has its own random number generator.
        output_handler->setupOutputForRun(current_run);
        int run_seed = current_run + current_sim_seed;
        const auto random_generator = std::make_unique<Randomizer>(run_seed);
        const auto site = createSites(current_run, random_generator.get(), analyser.get(), input_dir);
        site->setStoppingCondition(parameters_.stopping_criteria);
        SimulationTime time{parameters_.time_stepping, parameters_.max_time};

        // Start simulation for-loop over all timesteps for one run
        for (time.updateTimestep(0); !time.endReached(); ++time) {
            // Output of the state at the current time of the run as xml files and visualize
            output_handler->outputCurrentConfiguration(*site, time, current_run);
            visualizer->visualizeCurrentConfiguration(*site, time, current_run);

            // All interactions and dynamics of the hABM for all cells is performed
            site->doAgentDynamics(random_generator.get(), time);
            // Check if timestep size is changed
            site->updateTimeStepSize(time);
            if (site->checkForStopping(time)) {
                break;
            }
            if (time.checkForNumberOfExecutions(100, true)) {
                DEBUG_STDOUT("Run: " << current_run << ", Time: " << time.getCurrentTime());
            }
        }
        // Outputs for last timesteps
        output_handler->outputCurrentConfiguration(*site, time, current_run, run_seed, true,parameters_.cmd_input_args);
        visualizer->visualizeCurrentConfiguration(*site, time, current_run, true);
    }

    // Write outputs
    const auto end = timer::now();
    analyser->outputAllMeasurements();
    output_handler->concludeSimulation(runs,
                                       seed,
                                       static_cast<double>(std::chrono::duration_cast<std::chrono::seconds>(
                                               end - start).count()));
}

std::unique_ptr<Site> Simulator::createSites(int run, Randomizer *random_generator,
                                             const Analyser *analyser,
                                             const std::string &input_dir) const {
    if (parameters_.site_parameters->type == "SphereSite") {
        return std::make_unique<SphereSite>(parameters_.site_parameters.get(),
                                            parameters_.time_stepping,
                                            parameters_.dimensions,
                                            input_dir,
                                            random_generator,
                                            analyser->generateMeasurement(std::to_string(run)));
    }
    if (parameters_.site_parameters->type == "AlveoleSite") {
        return std::make_unique<AlveoleSite>(parameters_.site_parameters.get(),
                                             parameters_.time_stepping,
                                             parameters_.dimensions,
                                             input_dir,
                                             random_generator,
                                             analyser->generateMeasurement(std::to_string(run)));
    }
    ERROR_STDERR("Unknown site with type name: " << parameters_.site_parameters->type);
    return std::unique_ptr<Site>();
}

void Simulator::handleCmdInputArgs(const std::unordered_map<std::string, std::string> &cmd_input_args) {
    // Parameters to be screened or from cmd input
    // Here: dc, sAEC, nOfM and nOfCon that are specified in parameter_screening option in config_*.json
    for (const auto&[key, value] : cmd_input_args) {
        if ("dc" == key) {
            parameters_.site_parameters->particle_manager_parameters.diffusion_constant = std::stod(value);
            SYSTEM_STDOUT("Set parameter: " << key << " = " << value);
            updateTimestepForDC(parameters_.site_parameters->particle_manager_parameters.diffusion_constant);
        }
        if ("sAEC" == key) {
            parameters_.site_parameters->particle_manager_parameters.molecule_secretion_per_cell = std::stod(value);
            SYSTEM_STDOUT("Set parameter: " << key << " = " << value);
        }
    }
    for (const auto &agent : parameters_.site_parameters->agent_manager_parameters.agents) {
        if ("Macrophage" == agent->type) {
            for (const auto&[key, value] : cmd_input_args) {
                if ("nOfM" == key) {
                    agent->number = std::stod(value);
                    SYSTEM_STDOUT("Set parameter: " << key << " = " << value);
                }
            }
        }
        if ("AspergillusFumigatus" == agent->type) {
            auto *af_parameters = static_cast<abm::util::SimulationParameters::AgentParameters *>(agent.get());
            for (const auto&[key, value] : cmd_input_args) {
                if ("nOfCon" == key) {
                    af_parameters->number = std::stoi(value);
                    SYSTEM_STDOUT("Set parameter: " << key << " = " << value);
                }
            }
        }
    }
}

void Simulator::updateTimestepForDC(double dc) {
    // Stable timesteps for corresponding diffusion coefficients dc based on previous stability analysis
    if (abm::util::approxEqual(dc, 60))
        parameters_.time_stepping = 0.05;
    else if (abm::util::approxEqual(dc, 600))
        parameters_.time_stepping = 0.005;
    else if (abm::util::approxEqual(dc, 6000))
        parameters_.time_stepping = 0.0005;
    else {
        parameters_.time_stepping = 2 / dc;
    }
    SYSTEM_STDOUT("Updated timestep to " << parameters_.time_stepping);
}
