//  Copyright by Christoph Saffer, Paul Rudolph, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#include <boost/filesystem.hpp>
#include <iostream>
#include <map>

#include "io/InputConfiguration.h"
#include "simulation/simulator.h"
#include "utils/io_util.h"
#include "utils/macros.h"
#include "utils/misc_util.h"
#include <omp.h>

int main(int argc, char** argv) {

    // Default location of config file (if no parameter is specified)
    boost::filesystem::path config_xml("../../config.json");
    if (argc > 1 && !(std::istringstream(argv[1]) >> config_xml)) {
        ERROR_STDERR("usage: " << argv[0] << " <config.json>");
        return 1;
    }
    if (!(boost::filesystem::exists(config_xml))) {
        ERROR_STDERR("Configuration File in " << config_xml.string() << " does not exist!");
        ERROR_STDERR("usage: " << argv[0] << " <config.json>");
        return 2;
    }

    // Change root directory for simulation to configuration location
    chdir(&config_xml.parent_path().c_str()[0]);

    // Read cmd inputs and screening parameters
    const auto parameters = abm::util::getMainConfigParameters(config_xml.filename().string());
    std::unordered_map<std::string, std::string> input_args = abm::util::handleCmdInputs(argc, argv);

    // Initialize parallelization
#if defined(_OPENMP)
    omp_set_num_threads(parameters.number_of_threads);
    DEBUG_STDOUT("OpenMP activated with " << parameters.number_of_threads << " Thread(s).");
#endif

    // Start simulation runs
    if (parameters.screening_parameters.empty()) {
        const auto simulator = std::make_unique<const Simulator>(parameters.config_path, input_args);
        simulator->executeRuns(parameters.runs, parameters.system_seed, parameters.output_dir, parameters.input_dir);
    } else {
        // Screening over all parameter combinations specified as sets in the configuration file <config.json>
        // For screening, the cartesian product of all the single parameters sets is generated
        const auto [parameter_names, value_combinations] = abm::util::calculateCartesianProd(parameters.screening_parameters);
        for (int sim = 0; sim < value_combinations.size(); ++sim) {
            std::stringstream sim_para{};
            for (int i = 0; i < parameter_names.size(); ++i) {
                sim_para << parameter_names[i] << value_combinations[sim][i] << "_";
                input_args[parameter_names[i]] = value_combinations[sim][i];
            }
            SYSTEM_STDOUT("Start " << parameters.runs << " runs of simulation " << sim + 1 << "/" << value_combinations.size());
            const auto simulator = std::make_unique<const Simulator>(parameters.config_path, input_args);
            simulator->executeRuns(parameters.runs, parameters.system_seed, parameters.output_dir, parameters.input_dir, sim, sim_para.str());
        }
    }

    return 0;
}
