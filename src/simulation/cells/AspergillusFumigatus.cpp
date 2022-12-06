//  Copyright by Christoph Saffer, Paul Rudolph, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#include "AspergillusFumigatus.h"
#include "simulation/Site.h"
#include "simulation/morphology/SphericalMorphology.h"
#include "utils/macros.h"

std::string AspergillusFumigatus::getTypeName() {
    return "AspergillusFumigatus";
}

void AspergillusFumigatus::handleInteractionEvent(InteractionEvent *ievent) {
    if (ievent->getNextState() == "Lysis") {
        auto aspState = getCellStateByName("Death");
        INFO_STDOUT("Death of AspergillusFumigatus");
        if (aspState != 0) {
            setState(aspState);
        }
    }

    if (ievent->getNextState() == "Phagocytose") {
        auto swelling = getCellStateByName("AspergillusSwelling");
    }
}

void AspergillusFumigatus::move(double timestep, double current_time) {
}

void AspergillusFumigatus::doMorphologicalChanges(double timestep, double current_time) {
}

void AspergillusFumigatus::setup(double time_delta,
                                 double current_time,
                                 abm::util::SimulationParameters::AgentParameters *parameters) {
    Cell::setup(time_delta, current_time, parameters);
}
