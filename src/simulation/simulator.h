//  Copyright by Christoph Saffer, Paul Rudolph, Sandra Timme, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#ifndef SIMULATOR_SIMULATOR_H_
#define SIMULATOR_SIMULATOR_H_

#include "utils/io_util.h"

class Site;
class Analyser;
class Randomizer;

namespace abm::test { std::string test_simulation(const std::string &config); }
class Simulator {
public:
    /// Class for starting simulations
    Simulator(std::string config_path, const std::unordered_map<std::string, std::string> &cmd_input_args);
    ~Simulator();
    Simulator(const Simulator &) = delete;
    Simulator &operator=(const Simulator &) = delete;
    Simulator(Simulator &&) = delete;
    Simulator &operator=(Simulator &&) = delete;

    /**
     * Initializes all objects needed for a simulation and contains the main for-loop over all timesteps
     * @param runs Integer that contains the number of runs
     * @param seed Integer that contains the seed value
     * @param output_dir String that contains the output directory
     * @param input_dir String that contains the input directory
     * If parameter screening is activated:
     * @param sim Integer that contains the enumerated parameter configuration that is currently used
     * @param parameter_string String that contains the parameter configuration that is currently used
     */
    void executeRuns(int runs, int seed, const std::string &output_dir, const std::string &input_dir, int sim = 0,
                     const std::string &parameter_string = "") const;

    /*!
     * Creates the environment for the simulation
     * @param run Integer that contains the current run
     * @param random_generator Randomizer object
     * @param analyser Analyser object
     * @param input_dir String that contains input directory
     * @return Object of created Site (e.g. AlveoleSite)
     */
    std::unique_ptr<Site> createSites(int run,
                                      Randomizer *random_generator,
                                      const Analyser *analyser,
                                      const std::string &input_dir) const;

    /*!
     * Handles inputs from the command line OR from the option "parameter_screening" in the config.json file
     * @param cmd_input_args map of two Strings that contain the parameter and the assigned value
     */
    void handleCmdInputArgs(const std::unordered_map<std::string, std::string> &cmd_input_args);

    /*!
     * Adjusts the timestep according to the chosen diffusion coefficient dc
     * @param dc Double that contains diffusion coefficient dc
     */
    void updateTimestepForDC(double dc);

    /// Used for integration tests
    friend std::string abm::test::test_simulation(const std::string &config);

private:
    static int consumers;
    std::string config_path_{};
    abm::util::SimulationParameters parameters_{};
};

#endif // SIMULATOR_SIMULATOR_H_
