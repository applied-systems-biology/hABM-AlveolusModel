//  Copyright by Christoph Saffer, Paul Rudolph, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#ifndef VISUALISATIONMAIN_H
#define    VISUALISATIONMAIN_H

#include <string>
#include "utils/io_util.h"
#include "utils/time_util.h"


class Site;

class Visualizer {
public:
    /// Central class for visualizing a simulation
    Visualizer() = default;
    Visualizer(const std::string &config_path, const std::string &project_dir, int total_runs = 0);

    /*!
     * Visualize current configuration per timestep
     * @param site Site object of environment (e.g. AlveoleSite)
     * @param time SimulationTime object, i.e. contains current time and timestep
     * @param run Integer that contains current run
     * @param simulation_end Boolean that denotes if simulation end is reached
     */
    void visualizeCurrentConfiguration(const Site &site, const SimulationTime &time, int run,
                                       bool simulation_end = false) const;
private:
    void concludeRun() const;
    std::vector<std::string> visualization_path_{};
    mutable abm::util::VisualizerParameters parameters_{};
};

#endif    /* VISUALISATIONMAIN_H */

