//  Copyright by Christoph Saffer, Paul Rudolph, Sandra Timme, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#ifndef RANDOMWALK_H
#define	RANDOMWALK_H

#include "simulation/Movement.h"

#include "basic/Randomizer.h"
#include "basic/InversionSampler.h"

class RandomWalk : public Movement {
public:
  // Class for the normal random walk
    RandomWalk(Site *site, unsigned int spatial_dimensions);
    RandomWalk(Site *site, double speed, unsigned int spatial_dimensions);
    RandomWalk(Site *site, double speed, double speedStddev, unsigned int spatial_dimensions);

    Coordinate3D *move(double timestep, double dc);
    Coordinate3D *calculateDiffusiveMove(double, double dc);
    Coordinate3D *calculateRandomMove(double);
    double getSpeed() { return vector_length_per_timeunit_; };
    std::string getMovementName() final { return "RandomWalk"; }

private:
    double vector_length_per_timeunit_{};
    double speed_stddev_{};
};

#endif	/* RANDOMWALK_H */

