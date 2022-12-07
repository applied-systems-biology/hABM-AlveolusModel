//  Copyright by Christoph Saffer, Paul Rudolph, Sandra Timme, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#include <memory>

#include "simulation/RateFactory.h"
#include "io/InputConfiguration.h"
#include "simulation/Condition.h"
#include "simulation/rates/ConstantRate.h"
#include "simulation/rates/ConditionalRate.h"


std::unordered_map<std::string, std::unique_ptr<Rate>> RateFactory::rates_{};
bool RateFactory::initialized_ = false;

void RateFactory::initialize(const std::vector<RateParameter> &rates) {
    initialized_ = true;
    for (const auto &rate: rates) {
        if ("ConstantRate" == rate->type) {
            RateFactory::rates_.emplace(std::make_pair(rate->key, std::make_unique<ConstantRate>(rate->rate)));
        } else if ("ConditionalRate" == rate->type) {
            auto condition = std::make_unique<Condition>(rate->condition);
            RateFactory::rates_.emplace(
                    std::make_pair(rate->key, std::make_unique<ConditionalRate>(rate->rate, std::move(condition))));
        }
    }
}

void RateFactory::close() {
    rates_.clear();
}

std::unique_ptr<Rate> RateFactory::createRate(XMLNode rateNode) {
    return nullptr;
}

const Rate *RateFactory::getRate(const std::string &rate_key) {
    return RateFactory::rates_[rate_key].get();
}

const std::unordered_map<std::string, std::unique_ptr<Rate>> &RateFactory::getRates() {
    return RateFactory::rates_;
}
