//  Copyright by Christoph Saffer, Paul Rudolph, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#include "simulation/Movement.h"
#include "simulation/Site.h"


Movement::Movement(unsigned int spatial_dimensions) : spatial_dimensions_(spatial_dimensions) {
    current_move_ = std::make_unique<Coordinate3D>();
}

Coordinate3D *Movement::move(double ha, double diffusion_constant) {
    setCurrentTimestep(ha);
    return current_move_.get();
}

void Movement::setSite(Site *site) {
    site_ = site;
}

void Movement::setCurrentPosition(Coordinate3D *agent_pos) {
    current_pos_ = agent_pos;
}

void Movement::setCurrentTimestep(double timestep) {
    current_timestep_ = timestep;
}

void Movement::setCurrentMove(Coordinate3D *currMove) {
    *current_move_ = *currMove;
}

double Movement::getStartingTime() {
    return 0.0;
}

Coordinate3D *Movement::getCurrentMove() {
    return current_move_.get();
}

std::string Movement::getMovementName() {
    return std::string{};
}
