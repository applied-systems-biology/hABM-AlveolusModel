//  Copyright by Christoph Saffer, Paul Rudolph, Sandra Timme, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#include <cassert>

#include <boost/filesystem.hpp>

#include "Visualizer.h"
#include "utils/macros.h"
#include "utils/time_util.h"
#include "simulation/Site.h"
#include "visualisation/PovFile.h"
#include "visualisation/PovRayObject.h"

Visualizer::Visualizer(const std::string &config_path, const std::string &project_dir, int total_runs) {
    auto const visualizer_config = static_cast<boost::filesystem::path>(config_path).append(
            "visualisation-config.json").string();
    assert(!project_dir.empty());
    parameters_ = abm::util::getViualizerParameters(visualizer_config);

    // Initialize visualization folders for all runs
    if (parameters_.pov_active && total_runs > 0) {
        for (int i = 0; i < total_runs; i++) {
            std::string output_path = "povs-" + std::to_string(i + 1);
            visualization_path_.emplace_back(
                    static_cast<boost::filesystem::path>(project_dir).append(output_path).string());
            boost::filesystem::create_directories(visualization_path_[i]);
            DEBUG_STDOUT("Visualisation is activated for run id " << i + 1);
        }
    } else if (parameters_.pov_active && parameters_.run_id == 0) {
        ERROR_STDERR("Visualisation is activated but run id is 0");
    }
}

void Visualizer::visualizeCurrentConfiguration(const Site &site, const SimulationTime &time, int run,
                                               bool simulation_end) const {
    bool writeout = (0 == time % static_cast<int>(parameters_.output_interval));
    if (site.getParticleManager()->getDiffusionCoefficient() > 500) {
        // In ParticleManager::steadyStateReached - timestep changes during the simulation when a steady state is reached
        // Output interval must be adjusted to generate homogeneous output
        int frames_per_minute = parameters_.output_interval;
        double criteria_val = (frames_per_minute * time.getCurrentTime() -
                               round(frames_per_minute * time.getCurrentTime()));
        writeout = criteria_val < time.getCurrentDeltaT() / 2 && criteria_val > -time.getCurrentDeltaT() / 2;
    }

    // Generate visualization via PovFiles
    if (parameters_.pov_active && (writeout || simulation_end)) {
        const auto
                povDirTimestep = static_cast<boost::filesystem::path>(visualization_path_[run - 1]).append(
                abm::util::generateUniformString(time.getCurrentTime()));
        auto newPovFile = std::make_shared<PovFile>(povDirTimestep.string());
        const auto boundary = "";

        std::string globals;
        newPovFile->setGlobalPart(globals);
        newPovFile->setDimensions(parameters_.px_width, parameters_.px_height);
        newPovFile->setCameraPart(PovRayObject::getCamera(parameters_.camera_position, parameters_.camera_look_at,
                                                          parameters_.camera_angle));
        newPovFile->setBackgroundPart(PovRayObject::getBackground(ColorRGB(1, 1, 1)));
        newPovFile->setLightPart("");
        for (const auto &lightsource : parameters_.light_sources) {
            newPovFile->addLightPart(PovRayObject::getLightsource(lightsource));
        }
        newPovFile->setBorderPart(boundary);
        newPovFile->transcribeSite(site);
        newPovFile->transcribeParticles(site);
        newPovFile->transcribeAgents(site);
        newPovFile->doPovProcess(time.getCurrentTime());

        // Conclude simulation and generate video
        if (simulation_end) {
            concludeRun();
        }
    }
}

void Visualizer::concludeRun() const {
    using boost::filesystem::path;
    if (parameters_.output_video) {
        for (int i = 0; i < visualization_path_.size(); i++) {
            std::string w = parameters_.px_width, h = parameters_.px_height;
            std::ostringstream command;
            const auto log_pass = static_cast<path>(visualization_path_[i]).append("log.pass");
            command << "ffmpeg -hide_banner -loglevel panic -r 20.0 -pattern_type glob -i "
                    << static_cast<path>(visualization_path_[i]).append("*.png") << " -passlogfile " << log_pass
                    << " -an -vb 6k -s " << w << "x" << h << " -vbt 1M -pass 1 -threads 6  -vcodec png -f mp4 -y "
                    << static_cast<path>(visualization_path_[i]).append("action_movie.mp4") << " && "
                    << "ffmpeg -hide_banner -loglevel panic -r 20.0  -pattern_type glob -i "
                    << static_cast<path>(visualization_path_[i]).append("*.png") << " -passlogfile " << log_pass
                    << " -an -vb 6k -s " << w << "x" << h << "  -vbt 1M -pass 2 -threads 6  -vcodec png -f mp4 -y "
                    << static_cast<path>(visualization_path_[i]).append("action_movie.mp4");
            abm::util::executeShellCommand(command.str());
        }
    }
}

