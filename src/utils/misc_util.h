//  Copyright by Christoph Saffer, Paul Rudolph, Sandra Timme, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#ifndef ABM_UTILS_MISC_UTIL_H_
#define ABM_UTILS_MISC_UTIL_H_

#include <vector>
#include <string>
#include <memory>
#include <unordered_map>

#include "basic/Coordinate3D.h"
#include "basic/SphericCoordinate3D.h"

class Agent;
class Randomizer;
namespace abm::util {
    template<typename... Ts>
    struct overloaded : Ts ... {
        using Ts::operator()...;
    };
    template<typename... Ts> overloaded(Ts...) -> overloaded<Ts...>;
    template<char T = ';', typename... Args>
    std::string concatenate(Args &&... args);

    /*!
     * Generates hash from positions of outputs to compare simulation results
     * @param current_time Double for current time
     * @param agents vector of Agent that contains agents of current simulations
     * @return String that contains calculated hash value
     */
    std::string generateHashFromAgents(const double &current_time, const std::vector<std::shared_ptr<Agent>> &agents);

    /*!
     * Get file names from a directory
     * @param path String that contains directory
     * @param fileMask String that contains a file mask to only get files of one type (i.e. *.csv)
     * @return vector of String that contains all file names
     */
    std::vector<std::string> getFileNamesFromDirectory(const std::string &path, const std::string &fileMask = "");

    /*!
     * Generates cartesian product of n input sets
     * @param para Unordered map of Strings of parameters and corresponding values
     * @return pair of vector of String that contains all combinations of given parameters values
     */
    std::pair<std::vector<std::string>, std::vector<std::vector<std::string>>>
    calculateCartesianProd(const std::unordered_map<std::string, std::vector<std::string>> &para);

    /*!
     * Read coordinates from file
     * @param AMpos vector of Coordinate3D that contains written positions from file
     * @param inputString String that contains path to file
     */
    void read3DCoordinatesFromFile(std::vector<Coordinate3D> &AMpos, const std::string &inputString);

    /*!
     * read Lambda value from file for input of AM
     * @param input_string String that contains path to file
     * @return Double that contains Lambda value for input of AM
     */
    double readLambdaValueFromFile(const std::string &input_string);

    /*!
     * Checks if a folder exists
     * @return String that contains path to folder
     */
    bool folderExists(std::string);

    /*!
     * Checks if two values are approximately equal up to a certain accuracy
     * @param d1 Double that contains first number
     * @param d2 Double that contains second number
     * @param epsilon Double that contains accuracy
     * @return Bool if values are equal or not
     */
    bool approxEqual(double d1, double d2, double epsilon = 1e-8);

    /// Transforms SphericCoordinate3D to Coordinate3D
    Coordinate3D toCartesianCoordinates(const SphericCoordinate3D &coordinate);

    /// Transforms Coordinate3D to SphericCoordinate3D
    SphericCoordinate3D toSphericCoordinates(const Coordinate3D &coordinate);

    /// Rotates a Coordinate3D about an angle theta and phi
    Coordinate3D rotateAboutThetaAndPhi(const Coordinate3D &coordinate, double theta, double phi);

    /// Handles cmd inputs
    std::unordered_map<std::string, std::string> handleCmdInputs(int argc, char **argv);
}
#endif //ABM_UTILS_MISC_UTIL_H_
