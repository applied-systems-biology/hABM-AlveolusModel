//  Copyright by Christoph Saffer, Paul Rudolph, Sandra Timme, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#include "BoundaryCondition.h"
#include "simulation/Site.h"
#include "analyser/Analyser.h"


BoundaryCondition::BoundaryCondition(Site *site) {
    this->site = site;
}

void BoundaryCondition::handleBoundaryCross(Agent *agent, Coordinate3D *moveVec, double current_time) {
    return;
}

std::string BoundaryCondition::getTypeName() {
    return "BoundaryCondition";
}