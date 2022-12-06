//  Copyright by Christoph Saffer, Paul Rudolph, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#ifndef RATEFACTORY_H
#define RATEFACTORY_H

#include <boost/filesystem.hpp>

#include "io/XMLFile.h"
#include "utils/io_util.h"

class Rate;
class RateFactory {
    using RateParameter = std::unique_ptr<abm::util::InputParameters::DefaultRateParameters>;
public:
    /// Factory class for rate related calculations defining interaction rates. Rate definitions are given in the associated input-config.
    RateFactory() = delete;

    static void initialize(const std::vector<RateParameter> &rates);
    static void close();
    static std::unique_ptr<Rate> createRate(XMLNode rateNode);
    static const Rate *getRate(const std::string &rate_key);
    static const std::unordered_map<std::string, std::unique_ptr<Rate>> &getRates();
    static bool isInitialized() { return initialized_; };
private:
    static bool initialized_;
    static std::unordered_map<std::string, std::unique_ptr<Rate>> rates_;
};

#endif    /* RATEFACTORY_H */

