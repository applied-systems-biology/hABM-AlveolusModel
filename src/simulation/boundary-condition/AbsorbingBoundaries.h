//  Copyright by Christoph Saffer, Paul Rudolph, Sandra Timme, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#ifndef ABSORBINGBOUNDARIES_H
#define ABSORBINGBOUNDARIES_H

#include <memory>

#include "simulation/neighbourhood/BoundaryCondition.h"

class Analyser;

class AbsorbingBoundaries : public BoundaryCondition {
public:
  // Class for absorbing boundary conditions. This means that cells are removed from the system after they leave the environment. New cells are only inserted according to site-specific event.
    AbsorbingBoundaries(Site *site);

    void handleBoundaryCross(Agent *agent, Coordinate3D *moveVec, double current_time);
    std::string getTypeName();

};

#endif    /* ABSORBINGBOUNDARIES_H */

