//  Copyright by Christoph Saffer, Paul Rudolph, Sandra Timme, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#include "simulation/CellFactory.h"
#include "simulation/cells/Macrophage.h"
#include "simulation/cells/AspergillusFumigatus.h"
#include "simulation/Site.h"
#include "io/InputConfiguration.h"


std::map<double, XMLNode> CellFactory::donor_cells_{};
std::map<std::string, std::shared_ptr<abm::util::SimulationParameters::AgentParameters>> CellFactory::agent_configurations_{};

std::shared_ptr<Cell> CellFactory::createCell(const std::string &agenttype,
                                              std::unique_ptr<Coordinate3D> c,
                                              int id,
                                              Site *site,
                                              double time_delta,
                                              double current_time) {
    std::shared_ptr<Cell> agent{};
    auto site_tag = site->getIdentifier();
    if (agenttype == "Cell") {
        agent = std::make_shared<Cell>(std::move(c), id, site, time_delta, current_time);
    } else if (agenttype == "Macrophage") {
        agent = std::make_shared<Macrophage>(std::move(c), id, site, time_delta, current_time);
    } else if (agenttype == "AspergillusFumigatus") {
        agent = std::make_shared<AspergillusFumigatus>(std::move(c), id, site, time_delta, current_time);
    }
    agent->setup(time_delta, current_time, agent_configurations_[site_tag + agenttype].get());
    return agent;
}

void CellFactory::close() {
    agent_configurations_.clear();
}

void CellFactory::initialize(const std::unique_ptr<abm::util::SimulationParameters::SiteParameters> &site_parameters) {
    for (const auto &agent: site_parameters->agent_manager_parameters.agents) {
        agent_configurations_.emplace(std::make_pair(site_parameters->identifier + agent->type, agent));
    }
}

