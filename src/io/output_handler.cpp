//  Copyright by Christoph Saffer, Paul Rudolph, Sandra Timme, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#include <iostream>
#include <utility>
#include <boost/range/iterator_range.hpp>

#include "output_handler.h"
#include "utils/time_util.h"
#include "utils/io_util.h"
#include "io/XMLFile.h"
#include "utils/macros.h"
#include "simulation/Site.h"
#include "utils/misc_util.h"

using boost::filesystem::path;
using boost::filesystem::exists;

OutputHandler::OutputHandler(const std::string &config_path, std::string project_dir)
        : current_project_string_(static_cast<path>(project_dir).filename().string()),
          current_project_dir_(std::move(project_dir)) {

    parameters_ = abm::util::getOutputParameters(static_cast<path>(config_path).append("output-config.json").string());
    try {
        boost::filesystem::create_directories(current_project_dir_);
        //copy configuration files ending with "-config.xml"
        for (const auto &entry : boost::make_iterator_range(boost::filesystem::directory_iterator(config_path), {})) {
            if (const auto file_name = entry.path().filename().string(); file_name.find("config.json") !=
                                                                         std::string::npos) {
                boost::filesystem::copy(entry, static_cast<path>(current_project_dir_).append(file_name));
            }
        }

    } catch (const std::exception &e) {
        throw;
    }

}

void OutputHandler::concludeSimulation(int runs, int seed, double runtime) const {
    const std::string kRunsFile = static_cast<path>(current_project_dir_).parent_path().append("runs.xml").string();
    std::unique_ptr<XMLFile> xml_file;
    if (exists(kRunsFile)) {
        xml_file = std::make_unique<XMLFile>(kRunsFile, "Agent-Based-Framework");
    } else {
        xml_file = std::make_unique<XMLFile>("Agent-Based-Framework");
    }
    XMLNode run_node = xml_file->addChildToRootNode("Run");
    xml_file->addDataFieldToNode(run_node, "id", "discrete", "string", current_project_dir_);
    xml_file->addDataFieldToNode(run_node, "end date", "discrete", "date- and timestamp",
                                 abm::util::getCurrentLocalTimeAsString());
    xml_file->addDataFieldToNode(run_node, "configdir", "discrete", "string", current_project_dir_);
    xml_file->addDataFieldToNode(run_node, "tarfile", "discrete", "string",
                                 current_project_string_ + "-configs-*.tar.gz");
    xml_file->addDataFieldToNode(run_node, "number_of_runs", "discrete", "unsigned int", std::to_string(runs));
    xml_file->addDataFieldToNode(run_node, "runtime", "continuous", "double", std::to_string(runtime));
    xml_file->addDataFieldToNode(run_node, "randomizerseed", "discrete", "unsigned int", std::to_string(seed));

    std::ifstream buffStream{static_cast<path>(current_project_dir_).append("runs.csv").string()};
    std::vector<std::string> hashes;
    for (std::string line; std::getline(buffStream, line);) {
        hashes.emplace_back(line.substr(line.find_last_of(';') + 1, line.length()));
    }
    std::sort(hashes.begin(), hashes.end());
    std::ostringstream buff_string;
    std::copy(hashes.begin(), hashes.end(), std::ostream_iterator<std::string>(buff_string, " "));

    const auto hash_of_runs = std::hash<std::string>{}(buff_string.str());
    DEBUG_STDOUT("Hash: " << hash_of_runs);
    xml_file->addDataFieldToNode(run_node, "hash", "discrete", "unsigned long", std::to_string(hash_of_runs));
    xml_file->writeFile(kRunsFile);
}

void OutputHandler::concludeSimulationRun(int current_run, int seed, const std::string &hash,
                                          std::unordered_map<std::string, std::string> cmd_input_args) const {

    DEBUG_STDOUT("Run: " << current_run << ", Hash: " << hash);
    const auto current_run_string = current_project_string_ + "-configs-" + std::to_string(current_run);
    const auto current_run_tar = static_cast<path>(current_project_dir_).append(current_run_string + ".tar.gz");
    const auto current_runs_csv = static_cast<path>(current_project_dir_).append("runs.csv").string();
    const auto current_run_dir = static_cast<path>(current_project_dir_).append(current_run_string);
    std::ofstream csv_file(current_runs_csv, std::ofstream::out | std::ofstream::app);
    csv_file << current_run_string << ";" << seed << ';' << hash << '\n';
    csv_file.close();
    if (exists(current_run_dir)) {
        std::ostringstream cmd_compression;
        try {
            cmd_compression << "tar czf " << current_run_tar << " -C " << current_project_dir_ << " "
                            << current_run_string
                            << " --remove-files &";
            abm::util::executeShellCommand(cmd_compression.str());
        } catch (const std::exception &e) {
            throw;
        }
    }

    if (cmd_input_args.size() > 0) {
        if (current_run == 1) {
            const auto current_screening_csv = static_cast<boost::filesystem::path>(current_project_dir_).append(
                    "screening.csv").string();
            std::ofstream csv_file(current_screening_csv, std::ofstream::out | std::ofstream::app);
            csv_file << "seed,";
            for (const auto&[key, value] : cmd_input_args) {
                csv_file << key << ",";
            }

            csv_file << "\n" << seed << ",";
            for (const auto&[key, value] : cmd_input_args) {
                csv_file << value << ",";
            }
            csv_file.close();
        }
    }

}
void OutputHandler::setupOutputForRun(int run) const {
    if (parameters_.acitvated) {
        std::ostringstream run_directory;
        run_directory << current_project_string_ << "-configs-" << run;
#pragma omp critical
        {
            toc_dir_cache_[run] = static_cast<path>(current_project_dir_).append(run_directory.str()).append(
                    "toc").string();
            boost::filesystem::create_directories(toc_dir_cache_[run]);
        }
    }
}

void OutputHandler::outputCurrentConfiguration(const Site &site,
                                               const SimulationTime &time,
                                               int run,
                                               int seed,
                                               bool simulation_end,
                                               std::unordered_map<std::string, std::string> cmd_input_args) const {

    bool writeout = (0 == time % static_cast<int>(parameters_.output_interval));
    if (site.getParticleManager()->getDiffusionCoefficient() > 500) {
        // In ParticleManager::steadyStateReached - timestep changes during the simulation when a steady state is reached
        // Output interval must be adjusted to generate homogeneous output
        int frames_per_minute = parameters_.output_interval;
        double criteria_val = (frames_per_minute * time.getCurrentTime() -
                               round(frames_per_minute * time.getCurrentTime()));
        writeout = criteria_val < time.getCurrentDeltaT() / 2 && criteria_val > -time.getCurrentDeltaT() / 2;
    }
    // Generate xml output
    if (parameters_.acitvated && (writeout || simulation_end)) {
        const auto xml_tags = std::make_unique<XMLFile>();
        auto root_node = xml_tags->getRootNode();
        const auto name = abm::util::generateUniformString(time.getCurrentTime());
        xml_tags->addDataFieldToNode(root_node, "time", "discrete", "double",
                                     std::to_string(time.getCurrentTime()));
        site.includeSiteXMLTagToc(xml_tags.get());
        site.particle_manager_->includeParticleXMLTagToc(xml_tags.get());
        const auto csv_tags = site.agent_manager_->getAgentXMLTocTags(xml_tags.get(), parameters_.output_csv);
        xml_tags->writeFile(static_cast<path>(toc_dir_cache_.at(run)).append(name + ".xml").string());
        if (parameters_.output_csv) {
            std::ofstream csv_file(static_cast<path>(toc_dir_cache_.at(run)).append(name + ".csv").string());
            for (const auto &line : csv_tags) {
                csv_file << line << "\n";
            }
            csv_file.close();
        }
    }
    // Conclude simulation
    if (simulation_end && !current_project_dir_.empty()) {
        concludeSimulationRun(run, seed, abm::util::generateHashFromAgents(time.getCurrentTime(),
                                                                           site.agent_manager_->getAllAgents()), cmd_input_args);
    }
}