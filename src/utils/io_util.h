//  Copyright by Christoph Saffer, Paul Rudolph, Sandra Timme, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#ifndef ABM_UTILS_IO_UTIL_H_
#define ABM_UTILS_IO_UTIL_H_

#include <fstream>
#include <iostream>
#include <string>
#include <unordered_set>
#include <map>
#include <memory>
#include <unordered_map>
#include <array>
#include <type_traits>

#include <boost/filesystem.hpp>

#include "basic/Coordinate3D.h"

namespace abm::util {

    struct ConfigParameters {
        int runs{};
        int number_of_threads{};
        int system_seed{};
        std::string config_path{};
        std::string output_dir{};
        std::string input_dir{};
        std::unordered_map<std::string, std::vector<std::string>> screening_parameters{};
    };

    struct OutputParameters {
        int output_interval{};
        bool acitvated{};
        bool output_csv{};
    };

    struct SimulationParameters {
        struct MovementParameters {
            double diffusion_coefficient{};
            double persistence_time{};
            double mean{};
            double stddev{};
            std::string type{};
        };
        struct MorphologyParameters {
            double radius{};
            double stddev{};
            std::string type{};
            std::string color{};
        };
        struct CellStateParameters {
            std::string name{};
            std::vector<std::pair<std::string, std::string>> next_states;
        };

        struct BinomialDistribution {
            bool activated{};
            std::uint64_t n{};
            double p{};
        };
        struct MoleculeInteractionParameters {
            double secretion{};
            double uptake{};
            double k_blr{}; // ligand-receptor binding rate
            bool chemotaxis{};
            double k_r{}; //receptors recycling rate
            double k_i{}; //internalization rate
            double inside_conc_decay{}; // decay for conc inside the cell
        };

        struct AgentParameters {
            int initial_distribution{};
            float number{};
            std::string input_distribution_path{};
            std::string type{};
            BinomialDistribution binomial_distribution{};
            MorphologyParameters morphology_parameters{};
            MovementParameters movement_parameters{};
            MovementParameters passive_movement_parameters{};
            std::vector<CellStateParameters> states;
            std::vector<std::pair<std::string, MoleculeInteractionParameters>> molecule_interactions{};
        };

        struct ParticleManagerParameters {
            double diffusion_constant{};
            double molecule_secretion_per_cell{};
            bool draw_isolines{};
            std::string particle_delauney_input_file{};
        };

        struct MacrophageParameters : public AgentParameters {
            double k_r{};
            double k_i{};
            double k_blr{};
        };

        struct AgentManagerParameters {
            std::string site_identifier{};
            std::vector<std::shared_ptr<AgentParameters>> agents;
        };

        struct NHLParameters {
            std::string type{};
            int interaction_check_interval{};
            //for baloonlist
            double grid_constant{};
            double threshold{};
        };
        struct SiteParameters {
            bool passive_movement{};
            std::string identifier{};
            std::string type{};
            NHLParameters nhl_parameters{};
            ParticleManagerParameters particle_manager_parameters{};
            AgentManagerParameters agent_manager_parameters{};
        };
        struct AlveolusSiteParameter : public SiteParameters {
            int organism{};
            int objects_per_row{};
            int number_of_pok{};
            int number_of_aec2{};
            double site_radius{};
            double theta_lower_bound{};
            double r0_chemotaxis_source_aec1{};
            double r0_chemotaxis_source_aec2{};
            double lambda_input_rate{};
            double surfactant_thickness{};
            double thickness_of_border{};
            double radius_pores_of_kohn{};
            double radius_alv_epith_type_one{};
            double length_alv_epth_type_two{};
            Coordinate3D site_center{};
        };
        struct InteractionStateParameters {
            bool adhere{};
            double must_overhead{};
            std::string name{};
            std::string interaction_type{};
            std::vector<std::pair<std::string, std::string>> next_states;
        };

        struct InteractionParameters {
            //name and vector for all states the condition is restricted to, if the vector is empty all states of all cell are ok for the condition
            std::string name{};
            std::string type{};
            std::vector<std::pair<std::string, std::vector<std::string>>> cell_conditions;
            std::vector<InteractionStateParameters> states{};
        };

        int dimensions{};
        bool use_interactions{};
        double max_time{};
        double time_stepping{};
        std::vector<std::string> stopping_criteria{};
        std::string topic{};
        std::vector<std::unique_ptr<InteractionParameters>> interaction_parameters;
        std::unique_ptr<SiteParameters> site_parameters;
        std::unordered_map<std::string, std::string> cmd_input_args{};
    };

    struct AnalyserParameters {
        std::unordered_set<std::string> active_measurements{};
        std::vector<std::string> cell_state_count{};
    };

    struct VisualizerParameters {
        int run_id{};
        int total_runs{};
        int output_interval{};
        bool pov_active{};
        bool image_noise{};
        bool output_video{};
        double white_noise{};
        double camera_angle{};
        std::string px_width{};
        std::string px_height{};
        Coordinate3D camera_position{};
        Coordinate3D camera_look_at{};
        std::vector<Coordinate3D> light_sources;
    };

    struct InputParameters {
        struct DefaultRateParameters {
            double rate{};
            std::string type{};
            std::string key{};
            // for conditional rate
            std::string condition{};
        };
        struct DistributionParameters {
            std::string type{};
            std::string key{};
            std::string source_file{};
            int number_of_distributions{};
            double lambda{};
        };
        std::vector<std::unique_ptr<DefaultRateParameters>> rates{};
        std::vector<DistributionParameters> distributions{};
    };

    InputParameters getInputParameters(const std::string &input_config);
    VisualizerParameters getViualizerParameters(const std::string &visualizer_config);
    AnalyserParameters getAnalyserParameters(const std::string &analyser_config);
    ConfigParameters getMainConfigParameters(const std::string &config_path);
    OutputParameters getOutputParameters(const std::string &output_config);
    SimulationParameters getSimulationParameters(const std::string &simulator_config);

    void executeShellCommand(const std::string &command, bool suppress_output = true);

}
#endif//ABM_UTILS_IO_UTIL_H_
