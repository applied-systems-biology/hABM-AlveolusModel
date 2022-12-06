//  Copyright by Christoph Saffer, Paul Rudolph, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#include <cmath>

#include "analyser/InSituMeasurements.h"
#include "simulation/Site.h"
#include "simulation/site/AlveoleSite.h"
#include "simulation/neighbourhood/Collision.h"
#include "simulation/Interaction.h"
#include "utils/macros.h"

InSituMeasurements::InSituMeasurements(std::unordered_set<std::string> active_measurements, const std::string &id)
        : active_measurements_(std::move(active_measurements)) {

    for (const auto &active_:active_measurements_) {
        // Remove everything after "%" in the active measurement string
        auto active = active_;
        int char_pos = active_.find("%");
        if (char_pos > 0) {
            active = active_.substr(0, char_pos);
        }

        if ("alveole-statistics" == active) {
            pair_measurements_["alveole-statistics"] = std::make_unique<PairMeasurement>(id, "time", "nOfC", "nOfM");
        } else if ("agent-statistics" == active) {
            pair_measurements_["agent-statistics"] = std::make_unique<PairMeasurement>(id, "time", "agent", "agentid",
                                                                                       "state", "radius", "x", "y", "z",
                                                                                       "cellpart", "cellpart_id");
        } else if ("environment" == active) {
            pair_measurements_["environment"] = std::make_unique<PairMeasurement>(id, "x", "y", "z", "radius_or_length",
                                                                                  "type", "additional");
        }
    }
}

void InSituMeasurements::observeMeasurements(const SimulationTime &time) {
    using std::make_pair;
    const auto &current_time = time.getCurrentTime();

    for (const auto &active_:active_measurements_) {

        // Code snippet to let your measurement only execute ever X-th simulation step
        // Name your measurement as "measurment%X" -> example: "agent-statistics%25" ... every 25th step
        // Default if no value is given: Every 10th step
        int every_x_step = 10;
        auto active = active_;
        int char_pos = active_.find("%");
        if (char_pos > 0) {
            every_x_step = stoi(active_.substr(char_pos + 1, active_.size()));
            active = active_.substr(0, char_pos);
        }
        bool do_measurement = (time.getCurrentTimeStep() % every_x_step) == 0;

        // Measurements ever X timestep
        if ("alveole-statistics" == active && do_measurement) {
            pair_measurements_["alveole-statistics"]->addValuePairs(current_time,
                                                                    site_->getAgentManager()->getAgentQuantity(
                                                                            "AspergillusFumigatus"),
                                                                    site_->getAgentManager()->getAgentQuantity(
                                                                            "Macrophage"));
        } else if ("agent-statistics" == active && do_measurement) {
            for (auto agent: site_->getAgentManager()->getAllAgents()) {
                auto cellparts = agent->getAgentProperties()->getMorphology()->getAllSpheresOfThis().front();
                pair_measurements_["agent-statistics"]->addValuePairs(current_time,
                                                                      agent->getTypeName(), agent->getId(),
                                                                      agent->getCurrentCellState()->getStateName(),
                                                                      cellparts->getRadius(),
                                                                      cellparts->getPosition().x,
                                                                      cellparts->getPosition().y,
                                                                      cellparts->getPosition().z,
                                                                      "Mothercell", cellparts->getDescription());

            }
        }

        // End of the simulation
        if (time.getMaxTime() - time.getCurrentTime() <= time.getCurrentDeltaT() or site_->stopSimulation) {
            if ("environment" == active) {
                double radAEC1 = site_->getFeatureValueByName("radiusAlvEpithTypeOne");
                double lengthAEC2 = site_->getFeatureValueByName("lengthAlvEpithTypeTwo");
                double radPOK = site_->getFeatureValueByName("radiusPoresOfKohn");
                double thetaBorder = site_->getLowerThetaBound();
                double radius_alv = site_->getRadius();
                double thickness = site_->getThicknessOfBorder();
                pair_measurements_["environment"]->addValuePairs(0.0, 0.0, 0.0, radius_alv, "InnerAlveolus",
                                                                 thetaBorder);
                pair_measurements_["environment"]->addValuePairs(0.0, 0.0, 0.0, radius_alv + thickness,
                                                                 "OuterAlveolus",
                                                                 thetaBorder);

                for (auto x: site_->getAECT1()) {
                    std::string type = "AEC1";
                    Coordinate3D pos = abm::util::toCartesianCoordinates(x);
                    pair_measurements_["environment"]->addValuePairs(pos.x, pos.y, pos.z, radAEC1, type,
                                                                     thetaBorder);
                }
                for (auto x: site_->getAECT2()) {
                    std::string type = "AEC2";
                    Coordinate3D pos = abm::util::toCartesianCoordinates(x);
                    pair_measurements_["environment"]->addValuePairs(pos.x, pos.y, pos.z, lengthAEC2, type,
                                                                     thetaBorder);
                }
                for (auto x: site_->getPOK()) {
                    std::string type = "PoK";
                    Coordinate3D pos = abm::util::toCartesianCoordinates(x);
                    pair_measurements_["environment"]->addValuePairs(pos.x, pos.y, pos.z, radPOK, type,
                                                                     thetaBorder);
                }
            }
        }
    }
}

void InSituMeasurements::writeToFiles(const std::string &output_dir) const {
    for (const auto&[name, measurement]: histogram_measurements_) {
        const auto file_name = static_cast<boost::filesystem::path>(output_dir).append(name + ".csv").string();
        if (!boost::filesystem::exists(file_name)) {
            std::ofstream file{file_name};

            const auto &keys = measurement->getKeys();
            file << "id" << HistogramMeasurement::delimeter;
            for (const auto &key: keys) {
                file << key << HistogramMeasurement::delimeter;
            }
            file << '\n';
            file << *measurement;
            file.close();
        } else {
            std::ofstream file{file_name, std::ios_base::app};
            file << *measurement;
            file.close();
        }
    }
    for (const auto&[name, measurement]: pair_measurements_) {
        const auto file_name = static_cast<boost::filesystem::path>(output_dir).append(name + ".csv").string();
        if (!boost::filesystem::exists(file_name)) {
            std::ofstream file{file_name};
            const auto &keys = measurement->getKeys();
            file << "id" << PairMeasurement::delimeter;
            for (const auto &key: keys) {
                file << key << PairMeasurement::delimeter;
            }
            file << '\n';
            file << *measurement;
            file.close();
        } else {
            std::ofstream file{file_name, std::ios_base::app};
            file << *measurement;
            file.close();
        }
    }
}