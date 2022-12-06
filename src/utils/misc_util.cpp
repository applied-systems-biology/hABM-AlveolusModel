//  Copyright by Christoph Saffer, Paul Rudolph, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#include <dirent.h>
#include <set>
#include <sstream>
#include <filesystem>
#include <boost/algorithm/string.hpp>
#include <bits/basic_string.h>

#include "utils/misc_util.h"
#include "basic/SphericCoordinate3D.h"
#include "simulation/Agent.h"
#include "utils/macros.h"

namespace abm::util {

    std::string generateHashFromAgents(const double &current_time, const std::vector<std::shared_ptr<Agent>> &agents) {
        std::set<std::string> sorted_content;
        std::transform(agents.begin(),
                       agents.end(),
                       std::inserter(sorted_content, sorted_content.begin()),
                       [current_time](const auto &agent) {
                           return std::to_string(agent->getId()) + "-" +
                                  std::to_string(agent->getSurface()->getAllSpheresOfThis().size()) +
                                  "-" + agent->getAgentCSVTagAoc(current_time);
                       });
//        for (auto x: sorted_content) {
//            std::cout << x << " ";
//        }
//        std::cout << "\n";
        std::ostringstream buff_string;
        std::copy(sorted_content.begin(), sorted_content.end(), std::ostream_iterator<std::string>(buff_string, " "));
        return std::to_string(std::hash<std::string>{}(buff_string.str()));
    }

    template<char T, typename... Args>
    std::string concatenate(Args &&... args) {
        std::ostringstream output;
        ((output << std::forward<Args>(args) << T), ...);
        return output.str();
    }

    std::vector<std::string> getFileNamesFromDirectory(const std::string &path, const std::string &fileMask) {
        std::vector<std::string> fileList;
        DIR *dirp = opendir(path.c_str());
        struct dirent *dp;
        while ((dp = readdir(dirp)) != nullptr) {
            if (static_cast<std::string>(dp->d_name).find(fileMask) != std::string::npos) {
                fileList.emplace_back(dp->d_name);
            }
        }
        std::sort(fileList.begin(),
                  fileList.end()); // apply some order on files to have equal vector on different machines
        closedir(dirp);
        return fileList;
    }

    bool folderExists(std::string folder){
        return std::filesystem::exists(folder);
    }

    void read3DCoordinatesFromFile(std::vector<Coordinate3D> &AMpos, const std::string &inputString) {
        std::ostringstream am_dist_path;
        am_dist_path << inputString;
        std::ifstream file(am_dist_path.str().c_str());
        std::string currentLine;
        getline(file, currentLine);
        while (getline(file, currentLine)) {
            std::vector<std::string> tokens;
            boost::algorithm::split(tokens, currentLine, boost::is_any_of(","));
            double x = atof(tokens.at(1).c_str());
            double y = atof(tokens.at(2).c_str());
            double z = atof(tokens.at(3).c_str());

            Coordinate3D receivedCoordinate{x, y, z};
            AMpos.emplace_back(receivedCoordinate);
        }
        file.close();
    }

    double readLambdaValueFromFile(const std::string &inputString) {
        std::ostringstream AM_dist_path;
        AM_dist_path << inputString;
        std::ifstream file(AM_dist_path.str().c_str());
        std::string currentLine;
        getline(file, currentLine);
        return std::stod(currentLine);
    }

    std::unordered_map<std::string, std::string> handleCmdInputs(int argc, char **argv) {

        std::unordered_map<std::string, std::string> inputArgs;

        if (argc > 2) {
            if (argc % 2 == 0) {
                for (int i = 2; i < argc; i = i + 2) {
                    std::ostringstream ss, ssVal;
                    ss << argv[i];
                    std::string curArg = ss.str();
                    std::string firstCharArg = curArg.substr(0, 1);
                    if (firstCharArg == "-") {
                        ssVal << argv[i + 1];
                        std::string curVal = ssVal.str();
                        std::string withoutFirstCharArg = curArg.substr(1, -1);
                        inputArgs[withoutFirstCharArg] = curVal;
                    } else {
                        ERROR_STDERR(
                                R"(usage hint: "./hABM <config.json> -id1 value1 -id2 value2 ... ")");
                        ERROR_STDERR("now stopping execution");
                        exit(1);
                    }
                }
            } else {
                ERROR_STDERR(
                        R"(usage hint: "./hABM <config.json> -id1 value1 -id2 value2 ... ")");
                ERROR_STDERR("now stopping execution");
                exit(1);
            }
        }

        return inputArgs;
    }

    std::pair<std::vector<std::string>, std::vector<std::vector<std::string>>>
    calculateCartesianProd(const std::unordered_map<std::string, std::vector<std::string>> &para) {
        std::vector<std::string> keys;
        std::vector<std::vector<std::string>> values;
        std::vector<int> its;
        std::vector<std::vector<std::string>> tuples;
        for (const auto &x: para) {
            keys.emplace_back(x.first);
            values.emplace_back(x.second);
            its.emplace_back(0);
        }

        bool keep_running = true;
        while (keep_running) {
            std::vector<std::string> valtuple;
            for (int i = 0; i < keys.size(); ++i) {
                valtuple.emplace_back(values[i][its[i]]);
            }
            tuples.emplace_back(valtuple);

            for (int j = 0; j < keys.size();) {
                if (its[j] < values[j].size() - 1) {
                    its[j] += 1;
                    j = keys.size();
                } else {
                    if (j == keys.size() - 1) {
                        keep_running = false;
                    }
                    its[j] = 0;
                    ++j;
                }
            }

        }
        return {keys, tuples};
    }

    bool approxEqual(double d1, double d2, double epsilon) {
        return std::abs(d1 - d2) < epsilon;
    }

    Coordinate3D toCartesianCoordinates(const SphericCoordinate3D &coordinate) {
        double x = coordinate.r * sin(coordinate.theta) * cos(coordinate.phi);
        double y = coordinate.r * sin(coordinate.theta) * sin(coordinate.phi);
        double z = coordinate.r * cos(coordinate.theta);

        return {x, y, z};
    }

    SphericCoordinate3D toSphericCoordinates(const Coordinate3D &coordinate) {
        double r, theta, phi;
        r = sqrt(coordinate.x * coordinate.x + coordinate.y * coordinate.y + coordinate.z * coordinate.z);
        if (r == 0) {
            return {0, 0, 0};
        }
        theta = acos(coordinate.z / r);
        if (coordinate.x > 0) {
            phi = atan(coordinate.y / coordinate.x);
        } else {
            if (coordinate.x == 0.0) {
                if (coordinate.y >= 0) {
                    phi = M_PI / 2;
                } else {
                    phi = -M_PI / 2;
                }
            } else {
                if (coordinate.y >= 0.0) {
                    phi = atan(coordinate.y / coordinate.x) + M_PI;
                } else {
                    phi = atan(coordinate.y / coordinate.x) - M_PI;
                }
            }
        }
        return {r, theta, phi};
    }

    Coordinate3D rotateAboutThetaAndPhi(const Coordinate3D &coordinate, double theta, double phi) {
        SphericCoordinate3D pos_sph_old = toSphericCoordinates(coordinate);
        return toCartesianCoordinates({pos_sph_old.r, pos_sph_old.theta + theta, pos_sph_old.phi + phi});
    }

}