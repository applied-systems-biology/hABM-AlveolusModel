//  Copyright by Christoph Saffer, Paul Rudolph, Sandra Timme, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#ifndef ABM_OUTPUT_OUTPUTHANDLER_H_
#define ABM_OUTPUT_OUTPUTHANDLER_H_

#include <string>
#include <unordered_map>

#include "utils/io_util.h"

class Site;
class SimulationTime;

class OutputHandler {
public:
    /// Class for generating output files
    OutputHandler() = default;
    OutputHandler(const std::string &config_path, std::string project_dir);

    /*!
     * Finalizes the entire simulation
     * @param runs Integer that contains number of runs
     * @param seed Integer that contains system seed
     * @param runtime Double that contains runtime
     */
    void concludeSimulation(int runs, int seed, double runtime) const;

    /*!
     * Setups the output for a single run
     * @param run Integer that contains run number
     */
    void setupOutputForRun(int run) const;

    /*!
     * Writes output for a current timestep
     * @param site Site object of environment (e.g. AlveoleSite)
     * @param time SimulationTime object, i.e. contains current time and timestep
     * @param run Integer that contains run numbe
     * @param seed Integer that contains system seed
     * @param simulation_end Boolean that denotes if simulation end is reached
     * @param cmd_input_args Object that contains cmd inputs
     */
    void outputCurrentConfiguration(const Site &site,
                                    const SimulationTime &time,
                                    int run,
                                    int seed = 0,
                                    bool simulation_end = false,
                                    std::unordered_map<std::string, std::string> cmd_input_args = {}) const;

private:
    void concludeSimulationRun(int current_run, int seed, const std::string &hash,
                               std::unordered_map<std::string, std::string> cmd_input_args) const;
    std::string current_project_string_{};
    std::string current_project_dir_{};
    abm::util::OutputParameters parameters_{};
    abm::util::ConfigParameters config_param_{};
    mutable std::unordered_map<int, std::string> toc_dir_cache_;
};
#endif //ABM_OUTPUT_OUTPUTHANDLER_H_
