//  Copyright by Christoph Saffer, Paul Rudolph, Sandra Timme, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#ifndef M_PI
#define M_PI    3.14159265358979323846f
#endif

#include <math.h>

#include "RandomWalk.h"
#include "simulation/Site.h"
#include "utils/macros.h"

RandomWalk::RandomWalk(Site *site, unsigned int spatial_dimensions) : Movement(spatial_dimensions) {
    this->site_ = site;
    vector_length_per_timeunit_ = 0.1;
    speed_stddev_ = 0;
}

RandomWalk::RandomWalk(Site *site, double vectorLengthPerTimeUnit, unsigned int spatial_dimensions) : Movement(
        spatial_dimensions) {

    this->vector_length_per_timeunit_ = vectorLengthPerTimeUnit;
    speed_stddev_ = 0;
    this->site_ = site;
}

RandomWalk::RandomWalk(Site *site, double vectorLengthPerTimeUnit, double speedStddev, unsigned int spatial_dimensions)
        : Movement(spatial_dimensions) {

    this->site_ = site;
    this->vector_length_per_timeunit_ = vectorLengthPerTimeUnit;
    this->speed_stddev_ = speedStddev;
}

Coordinate3D *RandomWalk::move(double timestep, double dc) {
    setCurrentTimestep(timestep);
    if (dc > 0) {
        calculateDiffusiveMove(timestep, dc);
    } else {
        calculateRandomMove(timestep);
    }
    return current_move_.get();
}

Coordinate3D *RandomWalk::calculateRandomMove(double timestep) {
    double sampledSpeed, length;
    sampledSpeed = site_->getRandomGenerator()->generateNormalDistributedValue(vector_length_per_timeunit_,
                                                                               speed_stddev_);
    length = sampledSpeed * timestep;
    Coordinate3D prevMove = *current_move_;
    *current_move_ = site_->generateRandomDirectionVector(*current_pos_, length);
    return current_move_.get();
}

Coordinate3D *RandomWalk::calculateDiffusiveMove(double timestep, double dc) {
    DEBUG_STDOUT("calculate diffusive walk");
    // method holds only for 2D spatial systems
    double sigma = sqrt(2 * dc * timestep);
    DEBUG_STDOUT(std::to_string(dc) + " " + std::to_string(sigma));
    double length = site_->getRandomGenerator()->generateReighlayDistributedValue(sigma);
    *current_move_ = site_->generateRandomDirectionVector(*current_pos_, length);

    return current_move_.get();
}