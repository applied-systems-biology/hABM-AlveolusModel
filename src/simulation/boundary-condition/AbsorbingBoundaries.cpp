//  Copyright by Christoph Saffer, Paul Rudolph, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#include "AbsorbingBoundaries.h"
#include "simulation/Cell.h"
#include "simulation/Site.h"
#include "analyser/InSituMeasurements.h"

AbsorbingBoundaries::AbsorbingBoundaries(Site *site)
        : BoundaryCondition(site) {}

void AbsorbingBoundaries::handleBoundaryCross(Agent *agent, Coordinate3D *moveVec, double current_time) {
    site->getAgentManager()->replaceAgent(site, agent, nullptr, current_time);
}

std::string AbsorbingBoundaries::getTypeName() {
    return "AbsorbingBoundaries";
}