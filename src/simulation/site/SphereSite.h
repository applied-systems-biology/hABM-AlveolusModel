//  Copyright by Christoph Saffer, Paul Rudolph, Sandra Timme, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#ifndef SPHERESITE_H
#define    SPHERESITE_H

#include <string>
#include <iostream>

#include "simulation/Site.h"
#include "basic/Coordinate3D.h"

class SphereSite : public Site {
public:
  // Class for spherical class environment
    SphereSite(double time_delta,
               unsigned int spatial_dimensions,
               Randomizer *random_generator,
               std::shared_ptr<InSituMeasurements> measurements) : Site(time_delta,
                                                                        spatial_dimensions,
                                                                        random_generator, std::move(measurements)) {};

    SphereSite(abm::util::SimulationParameters::SiteParameters *parameters,
               double time_delta,
               unsigned int spatial_dimensions,
               const std::string &input_dir,
               Randomizer *random_generator,
               std::shared_ptr<InSituMeasurements> measurements);

    void includeSiteXMLTagToc(XMLFile *xmlFile) const override;
    void handleBoundaryCross(Agent *, Coordinate3D *, double current_time) final;
    bool containsPosition(Coordinate3D) override;
    [[nodiscard]] double getRadius() const final { return radius; };
    [[nodiscard]] std::string getType() const override { return "SphereSite"; }

    Coordinate3D getRandomPosition() override;
    Coordinate3D getRandomBoundaryPoint() override;
    Coordinate3D getCenterPosition() final;
    Coordinate3D getLowerLimits() final;
    Coordinate3D getUpperLimits() final;
    Coordinate3D generateRandomDirectionVector(Coordinate3D position, double length) override;
    Coordinate3D generatePersistentDirectionVector(Coordinate3D position,
                                                   double length,
                                                   Coordinate3D prevVector,
                                                   double previousAlpha) override;
    Coordinate3D generateBackShiftOnContacting(SphereRepresentation *activeSphere,
                                               SphereRepresentation *passiveSphere,
                                               double mustOverhead) final;

protected:
    virtual void adjustAgents(double time_delta, double current_time);

    double radius;
    Coordinate3D centerOfSite;
};

#endif    /* SPHERESITE_H */

