//  Copyright by Christoph Saffer, Paul Rudolph, Sandra Timme, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#include "testConfigurations.h"

#include <memory>

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "external/doctest/doctest.h"
#include "simulation/Site.h"
#include "simulation/simulator.h"
#include "analyser/Analyser.h"

using boost::filesystem::path;
using boost::filesystem::exists;

std::string abm::test::test_simulation(const std::string &config) {
  const auto parameters = abm::util::getMainConfigParameters(config);
  auto simulator = std::make_unique<Simulator>(parameters.config_path, std::unordered_map<std::string, std::string> {});
  const auto run_seed = parameters.system_seed;
  const auto analyser = std::make_unique<Analyser>();
  const auto random_generator = std::make_unique<Randomizer>(run_seed);
  const auto sites = simulator->createSites(0, random_generator.get(), analyser.get(), parameters.input_dir);
  const auto &site = sites;
  SimulationTime time{simulator->parameters_.time_stepping, simulator->parameters_.max_time};
  for (time.updateTimestep(0); !time.endReached(); ++time) {    //inner loop: t -> t + dt
    site->doAgentDynamics(random_generator.get(), time);
    site->updateTimeStepSize(time);
    if (site->checkForStopping(time)) {
      break;
    }
  }
  return abm::util::generateHashFromAgents(time.getCurrentTime(), site->getAgentManager()->getAllAgents());
}

// Alveolus Model Test Human
TEST_CASE ("Check Alveolus Test Human") {
    path config("../../test/configurations/testAlveolusHuman/config.json");
    CHECK(exists(config) == true);
    const auto string_return = abm::test::test_simulation(config.string());
    CHECK(string_return == "16053041708414146768");
}

TEST_CASE ("Check Alveolus Mouse Test") {
    path config("../../test/configurations/testAlveolusMouse/config.json");
    CHECK(exists(config) == true);
    const auto string_return = abm::test::test_simulation(config.string());
    CHECK(string_return == "502083924311758751");
}
