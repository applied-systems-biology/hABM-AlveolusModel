//  Copyright by Christoph Saffer, Paul Rudolph, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#include "simulation/BalloonListNHLocator.h"
#include "SphereSite.h"
#include "basic/Randomizer.h"
#include "basic/SphericCoordinate3D.h"
#include "simulation/Interactions.h"
#include "simulation/ParticleManager.h"
#include "utils/macros.h"
#include "simulation/AgentManager.h"

SphereSite::SphereSite(abm::util::SimulationParameters::SiteParameters *parameters,
                       double time_delta,
                       unsigned int spatial_dimensions,
                       const std::string &input_dir,
                       Randomizer *random_generator,
                       std::shared_ptr<InSituMeasurements> measurements) : Site(time_delta,
                                                                                spatial_dimensions,
                                                                                random_generator,
                                                                                std::move(measurements)) {


    double gridConstant = parameters->nhl_parameters.grid_constant;
    neighbourhood_locator_ = std::make_unique<BalloonListNHLocator>(gridConstant, getLowerLimits(), getUpperLimits(), this);

    if (neighbourhood_locator_ != nullptr) {
        neighbourhood_locator_->setThresholdDistance(parameters->nhl_parameters.threshold);
        unsigned int icInterval = static_cast<unsigned int>(parameters->nhl_parameters.interaction_check_interval);
        neighbourhood_locator_->setInteractionCheckInterval(icInterval);
    }
    initializeAgents(parameters->agent_manager_parameters, input_dir, 0, time_delta);
    particle_manager_->initializeParticles(this, parameters->particle_manager_parameters, input_dir);
    passiveMovementOn = parameters->passive_movement;

}

Coordinate3D SphereSite::getRandomPosition() {
    double x, y, z, phi, theta;
    double u, subst;
    switch (dimensions) {
        case 2: //on/in the surface of the sphere site
            u = random_generator_->generateDouble();
            phi = random_generator_->generateDouble(M_PI * 2.0);

            subst = 2 * radius * sqrt(u * (1 - u));
            x = subst * cos(phi);
            y = subst * sin(phi);
            z = radius * (1 - 2 * u);
            break;
        case 3:
            do {
                x = random_generator_->generateDouble(-1.0 * radius, radius);
                y = random_generator_->generateDouble(-1.0 * radius, radius);
                z = random_generator_->generateDouble(-1.0 * radius, radius);
            } while (x * x + y * y + z * z >= radius * radius);
            break;
        default:
            do {
                x = random_generator_->generateDouble(-1.0 * radius, radius);
                y = random_generator_->generateDouble(-1.0 * radius, radius);
                z = random_generator_->generateDouble(-1.0 * radius, radius);
            } while (x * x + y * y + z * z >= radius * radius);
            break;
    }
    //offset of the given coordinate elements to the center position of this site
    x += centerOfSite.x;
    y += centerOfSite.y;
    z += centerOfSite.z;

    return Coordinate3D{x, y, z};
}

Coordinate3D SphereSite::getCenterPosition() {
    return centerOfSite;
}

void SphereSite::handleBoundaryCross(Agent *agent, Coordinate3D *moveVec, double current_time) {
    DEBUG_STDOUT("Handling crossing of " + agent->getTypeName() + " " +
                 std::to_string(agent->getId()) + " at position " + agent->getPosition().printCoordinates() +
                                                                  " at timepoint " +
                                                                  std::to_string(current_time));
    boundary_condition_->handleBoundaryCross(agent, moveVec, current_time);

}

Coordinate3D SphereSite::getLowerLimits() {
    return Coordinate3D{centerOfSite.x - 1.2 * radius, centerOfSite.y - 1.2 * radius, centerOfSite.z - 1.2 * radius};
}

Coordinate3D SphereSite::getUpperLimits() {
    return Coordinate3D{centerOfSite.x + 1.2 * radius, centerOfSite.y + 1.2 * radius, centerOfSite.z + 1.2 * radius};
}

Coordinate3D SphereSite::getRandomBoundaryPoint() {
    double x, y, z;
    double u, subst, phi;

    u = random_generator_->generateDouble();
    phi = random_generator_->generateDouble(M_PI * 2.0);

    subst = 2 * radius * sqrt(u * (1 - u));
    x = subst * cos(phi);
    y = subst * sin(phi);
    z = radius * (1 - 2 * u);


    //offset of the given coordinate elements to the center position of this site
    x += centerOfSite.x;
    y += centerOfSite.y;
    z += centerOfSite.z;

    return Coordinate3D{x, y, z};
}

bool SphereSite::containsPosition(Coordinate3D position) {
    bool contains;

    double xWoOffset, yWoOffset, zWoOffset;

    xWoOffset = position.x - centerOfSite.x;
    yWoOffset = position.y - centerOfSite.y;
    zWoOffset = position.z - centerOfSite.z;

    //workaround for initialisation of morphology elements to check if sphericalmorphology is inside the sphere --> radius+0.1 was chosen instead of radius
    contains = (xWoOffset * xWoOffset + yWoOffset * yWoOffset + zWoOffset * zWoOffset <=
                (radius + 0.1) * (radius + 0.1));

    return contains;
}

Coordinate3D SphereSite::generateRandomDirectionVector(Coordinate3D position, double length) {
    double x, y, z, subst, r, phi, theta, u, dtheta, alpha, vix, viy, viz;
    SphericCoordinate3D currentPos = abm::util::toSphericCoordinates(position - centerOfSite);

    switch (dimensions) {
        case 2:

            dtheta = length / radius; //in rad
            alpha = random_generator_->generateDouble(M_PI * 2.0); //direction of the vector

            vix = radius * sin(dtheta) * cos(alpha);
            viy = radius * sin(dtheta) * sin(alpha);
            viz = radius * (cos(dtheta) - 1);

            phi = currentPos.phi;
            theta = currentPos.theta;

            x = cos(phi) * (cos(theta) * vix + sin(theta) * viz) - sin(phi) * viy;
            y = sin(phi) * (cos(theta) * vix + sin(theta) * viz) + cos(phi) * viy;
            z = cos(theta) * viz - sin(theta) * vix;

            break;

        case 3:
            //get a sin() sampled value in [0,PI] for theta via a uniform value of u
            u = random_generator_->generateDouble();//sampler->sample();

            phi = random_generator_->generateDouble(M_PI * 2.0);
            r = length;

            subst = 2 * r * sqrt(u * (1 - u));
            x = subst * cos(phi);
            y = subst * sin(phi);
            z = r * (1 - 2 * u);

            break;

        default:
            //get a sin() sampled value in [0,PI] for theta via a uniform value of u
            u = random_generator_->generateDouble();//sampler->sample();

            phi = random_generator_->generateDouble(M_PI * 2.0);
            r = length;

            subst = 2 * r * sqrt(u * (1 - u));
            x = subst * cos(phi);
            y = subst * sin(phi);
            z = r * (1 - 2 * u);

            break;

    }

    return Coordinate3D{x, y, z};
}

void SphereSite::adjustAgents(double time_delta, double current_time) {

    //repositioning to be on the surface again
    auto &allAgents = agent_manager_->getAllAgents();
    auto it = allAgents.begin();

    Coordinate3D origCoord, newCoord;
    SphericCoordinate3D sphereCoord;

    while (it != allAgents.end()) {
        origCoord = (*it)->getPosition();
        sphereCoord = abm::util::toSphericCoordinates(origCoord - centerOfSite);
        sphereCoord.r = radius;
        newCoord = abm::util::toCartesianCoordinates(sphereCoord);
        (*it)->setPosition(newCoord);

        it++;
    }

    //resolve conflicts
    it = allAgents.begin();

    while (it != allAgents.end()) {
        (*it)->getAgentProperties()->getInteractions()->avoidNewInteractions(time_delta, current_time);
        it++;
    }

}

void SphereSite::includeSiteXMLTagToc(XMLFile *xmlTags) const {
    std::ostringstream ssRadius;
    ssRadius << radius;
    XMLNode siteNode = xmlTags->addChildToRootNode("Site");

    XMLNode sphereSiteNode = xmlTags->addChildToNode(siteNode, "SphereSite");
    xmlTags->addDataFieldToNode(sphereSiteNode, centerOfSite);
    xmlTags->addDataFieldToNode(sphereSiteNode, "radius", "discrete", "double", ssRadius.str());

}

Coordinate3D SphereSite::generateBackShiftOnContacting(SphereRepresentation *activeSphere,
                                                       SphereRepresentation *passiveSphere,
                                                       double mustOverhead) {
    Coordinate3D normalVecBetwCells = activeSphere->getEffectiveConnection(passiveSphere);

    /*
     * push the cell away to the contact position in spheric systems
     */
    Coordinate3D r2Pos = passiveSphere->getPosition();
    Coordinate3D r2Vec = Coordinate3D{r2Pos.x, r2Pos.y, r2Pos.z};
    double rx, ry, rz;
    rx = r2Pos.x;
    ry = r2Pos.y;
    rz = r2Pos.z;

    //calculate normal vector of the surface where the translation will happen
    Coordinate3D n = Coordinate3D();
    n = r2Vec;
    n = n.crossProduct(normalVecBetwCells);
    n.setMagnitude(1.0);// *= 1.0/n.getMagnitude();
    double nx, ny, nz;
    nx = n.x;
    ny = n.y;
    nz = n.z;

    //the desired distance between the spheres after handling
    double d = passiveSphere->getRadius() + activeSphere->getRadius() - mustOverhead;

    // solution one and two, analytical result from mathematica
    double x1, y1, z1, x2, y2, z2;
    double diskriminante = -d * d * (nz * ry - ny * rz) * (nz * ry - ny * rz)
                           * (d * d * (nx * nx + ny * ny + nz * nz) - 4 * (ny * ny + nz * nz) * rx * rx +
                              8 * nx * ny * rx * ry
                              - 4 * (nx * nx + nz * nz) * ry * ry + 8 * nz * (nx * rx + ny * ry) * rz -
                              4 * (nx * nx + ny * ny) * rz * rz);
    double nenner = 2
                    * (nz * nz * (rx * rx + ry * ry) - 2 * nx * nz * rx * rz - 2 * ny * ry * (nx * rx + nz * rz) +
                       ny * ny * (rx * rx + rz * rz)
                       + nx * nx * (ry * ry + rz * rz));
    double xPart = d * d * (-ny * ny * rx + nx * ny * ry + nz * (-nz * rx + nx * rz));
    double yPart = d * d * (-nz * ry + ny * rz) * (-nx * ny * rx + nx * nx * ry + nz * (nz * ry - ny * rz));
    double zPart = d * d * (nz * ry - ny * rz) * (nz * (nx * rx + ny * ry) - (nx * nx + ny * ny) * rz);

    x1 = (-sqrt(diskriminante) + xPart) / (nenner);
    y1 = ((nz * rx - nx * rz) * sqrt(diskriminante) + yPart) / ((nz * ry - ny * rz) * nenner);
    z1 = ((-ny * rx + nx * ry) * sqrt(diskriminante) + zPart) / ((nz * ry - ny * rz) * nenner);;
    Coordinate3D dis1 = Coordinate3D{x1, y1, z1};

    x2 = (sqrt(diskriminante) + xPart) / (nenner);
    y2 = ((-nz * rx + nx * rz) * sqrt(diskriminante) + yPart) / ((nz * ry - ny * rz) * nenner);
    z2 = ((ny * rx - nx * ry) * sqrt(diskriminante) + zPart) / ((nz * ry - ny * rz) * nenner);;
    Coordinate3D dis2 = Coordinate3D{x2, y2, z2};

    Coordinate3D newPos1 = r2Vec + dis1;
    Coordinate3D newPos11 = newPos1;

    Coordinate3D newPos2 = r2Vec + dis2;
    Coordinate3D newPos22 = newPos2;

    double distance1 = newPos11.calculateEuclidianDistance(activeSphere->getPosition());
    double distance2 = newPos22.calculateEuclidianDistance(activeSphere->getPosition());
    Coordinate3D backShift = Coordinate3D();
    if (distance1 < distance2) {
        backShift = newPos11 - activeSphere->getPosition();
    } else {
        backShift = (newPos22) - activeSphere->getPosition();
    }

    return backShift;
}
Coordinate3D SphereSite::generatePersistentDirectionVector(Coordinate3D position,
                                                           double length,
                                                           Coordinate3D prevVector,
                                                           double previousAlpha) {
    return {};
}