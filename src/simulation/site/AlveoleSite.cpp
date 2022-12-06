//  Copyright by Christoph Saffer, Paul Rudolph, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#include <iterator>
#include <boost/algorithm/string.hpp>

#include "AlveoleSite.h"
#include "utils/macros.h"
#include "analyser/InSituMeasurements.h"
#include "simulation/BalloonListNHLocator.h"
#include "simulation/CellFactory.h"
#include "simulation/ParticleManager.h"
#include "simulation/AgentManager.h"

AlveoleSite::AlveoleSite(abm::util::SimulationParameters::SiteParameters *parameters,
                         double time_delta,
                         unsigned int spatial_dimensions,
                         const std::string &input_dir,
                         Randomizer *random_generator,
                         std::shared_ptr<InSituMeasurements> measurements) : SphereSite(time_delta,
                                                                                        spatial_dimensions,
                                                                                        random_generator,
                                                                                        std::move(measurements)) {

    // All parameters are set for the human or murine AlveoleSite
    identifier_ = parameters->identifier;
    auto alveolus_parameters = static_cast<abm::util::SimulationParameters::AlveolusSiteParameter *>(parameters);
    setBoundaryCondition();
    centerOfSite = alveolus_parameters->site_center;
    radius = alveolus_parameters->site_radius;
    thetaLowerBound = alveolus_parameters->theta_lower_bound;
    siteRadiusMin = radius;
    r0AEC1 = alveolus_parameters->r0_chemotaxis_source_aec1;
    r0AEC2 = alveolus_parameters->r0_chemotaxis_source_aec2;
    inputRate = alveolus_parameters->lambda_input_rate;
    surfactantThickness = alveolus_parameters->surfactant_thickness;
    thicknessOfBorder = alveolus_parameters->thickness_of_border;
    radiusPoresOfKohn = alveolus_parameters->radius_pores_of_kohn;
    radiusAlvEpithTypeOne = alveolus_parameters->radius_alv_epith_type_one;
    lengthAlvEpithTypeTwo = alveolus_parameters->length_alv_epth_type_two;
    noOfPoK = alveolus_parameters->number_of_pok;
    noOfAEC2 = alveolus_parameters->number_of_aec2;
    opR = alveolus_parameters->objects_per_row;
    organism = alveolus_parameters->organism;

    // Build alveolus
    includeAlveolarEpitheliumType1();
    if (organism == 1) {
        includePoresOfKohnHuman();
        includeAlveolarEpitheliumType2Human();
    } else {
        includePoresOfKohnMouse();
        includeAlveolarEpitheliumType2Mouse();
    }

    // Initialize neighbourhood locator
    double gridConstant = parameters->nhl_parameters.grid_constant;
    neighbourhood_locator_ = std::make_unique<BalloonListNHLocator>(gridConstant, getLowerLimits(), getUpperLimits(), this);
    if (neighbourhood_locator_ != nullptr) {
        neighbourhood_locator_->setThresholdDistance(parameters->nhl_parameters.threshold);
        unsigned int icInterval = static_cast<unsigned int>(parameters->nhl_parameters.interaction_check_interval);
        neighbourhood_locator_->setInteractionCheckInterval(icInterval);
    }

    // Insert agents into alveolus
    initializeAgents(parameters->agent_manager_parameters, input_dir, 0, time_delta);
    agent_manager_->setInitConQuantity();
    adjustAgents(time_delta, 0);//

    // Initialize particles and basic variables
    this->getParticleManager()->setAECCells(alvEpithTypeOne, alvEpithTypeTwo);
    particle_manager_->initializeParticles(this, parameters->particle_manager_parameters, input_dir);
    state_ = 1;
    boundary_input_vector_ = Coordinate3D();
    passiveMovementOn = parameters->passive_movement;
}

bool AlveoleSite::overPOK(Coordinate3D position) {
    SphericCoordinate3D sPos = abm::util::toSphericCoordinates(position - centerOfSite);
    bool insideMainSite = (sPos.r <= radius + thicknessOfBorder / 2 && sPos.r >= radius - thicknessOfBorder / 2)
                          && (sPos.theta >= thetaLowerBound);
    bool insidePoreOfKohn = false;
    if (insideMainSite) {
        double distance;
        for (const auto &pore : poresOfKohn) {
            Coordinate3D poreOfKohnPos = abm::util::toCartesianCoordinates(pore);
            distance = position.calculateEuclidianDistance(poreOfKohnPos);
            if (distance < radiusPoresOfKohn) {
                insidePoreOfKohn = true;
                break;
            }
        }
    }
    return insidePoreOfKohn;
}

bool AlveoleSite::containsPosition(Coordinate3D position) {
    SphericCoordinate3D sPos = abm::util::toSphericCoordinates(position - centerOfSite);
    bool insideMainSite = (sPos.r <= radius + thicknessOfBorder / 2 && sPos.r >= radius - thicknessOfBorder / 2)
                          && (sPos.theta >= thetaLowerBound);
    bool insidePoreOfKohn = false;
    if (insideMainSite) {
        double distance;
        for (const auto &pore : poresOfKohn) {
            Coordinate3D poreOfKohnPos = abm::util::toCartesianCoordinates(pore);
            distance = position.calculateEuclidianDistance(poreOfKohnPos);
            if (distance < radiusPoresOfKohn) {
                insidePoreOfKohn = true;
                break;
            }
        }
    }
    bool inside = insideMainSite && !insidePoreOfKohn;
    return inside;
}

void AlveoleSite::includeSiteXMLTagToc(XMLFile *xmlTags) const {
    std::ostringstream ssRadius, ssThetaLowerBound, ssThickness;
    ssRadius << radius;
    ssThetaLowerBound << thetaLowerBound;
    XMLNode siteNode = xmlTags->addChildToRootNode("Site");
    ssThickness << thicknessOfBorder;
    XMLNode sphereSiteNode = xmlTags->addChildToNode(siteNode, "AlveoleSite");
    xmlTags->addDataFieldToNode(sphereSiteNode, centerOfSite);
    xmlTags->addDataFieldToNode(sphereSiteNode, "radius", "continuous", "double", ssRadius.str());
    xmlTags->addDataFieldToNode(sphereSiteNode, "theta-lower-bound", "continuous", "double", ssThetaLowerBound.str());
    xmlTags->addDataFieldToNode(sphereSiteNode, "thickness", "continuous", "double", ssThickness.str());
    XMLNode currentEpithelialCellNode;

    for (auto it = alvEpithTypeOne.begin(); it != alvEpithTypeOne.end(); it++) {
        currentEpithelialCellNode = xmlTags->addChildToNode(sphereSiteNode, "AET1");
        Coordinate3D curPos = abm::util::toCartesianCoordinates(*it);
        xmlTags->addDataFieldToNode(currentEpithelialCellNode, curPos);
        std::ostringstream ssRadiusEpT;
        ssRadiusEpT << radiusAlvEpithTypeOne;
        xmlTags->addDataFieldToNode(currentEpithelialCellNode, "radius", "discrete", "double", ssRadiusEpT.str());
        XMLNode currentEpithelialNeigbourCellNode;
        //add all neighbouring type I cells
        for (auto it2 = alvEpithTypeOne.begin(); it2 != alvEpithTypeOne.end(); it2++) {
            Coordinate3D dConnect = Coordinate3D(
                    abm::util::toCartesianCoordinates(*it2) - abm::util::toCartesianCoordinates(*it));
            double distance = dConnect.getMagnitude();
            if (distance < 2 * radiusAlvEpithTypeOne && distance > 0.001) {
                currentEpithelialNeigbourCellNode = xmlTags->addChildToNode(currentEpithelialCellNode, "AETNC1");
                SphericCoordinate3D dConnectSphere = abm::util::toSphericCoordinates(dConnect);
                xmlTags->addDataFieldToNode(currentEpithelialNeigbourCellNode, dConnectSphere);
            }
        }
        //add all neighbouring type II cells
        for (auto it2 = alvEpithTypeTwo.begin(); it2 != alvEpithTypeTwo.end(); it2++) {
            Coordinate3D dConnect = Coordinate3D(
                    abm::util::toCartesianCoordinates(*it2) - abm::util::toCartesianCoordinates(*it));
            double distance = dConnect.getMagnitude();
            if (distance < sqrt(3.0) * lengthAlvEpithTypeTwo + radiusAlvEpithTypeOne) {
                currentEpithelialNeigbourCellNode = xmlTags->addChildToNode(currentEpithelialCellNode, "AETNC2");
                Coordinate3D neighbourBox = abm::util::toCartesianCoordinates(*it2);
                std::ostringstream ssLengthEpT;
                ssLengthEpT << lengthAlvEpithTypeTwo;
                xmlTags->addDataFieldToNode(currentEpithelialNeigbourCellNode, "length", "discrete", "double",
                                            ssLengthEpT.str());
                xmlTags->addDataFieldToNode(currentEpithelialNeigbourCellNode, neighbourBox);
            }
        }
        //add all neighbouring Pores of Kohn
        for (auto it2 = poresOfKohn.begin(); it2 != poresOfKohn.end(); it2++) {
            Coordinate3D dConnect = Coordinate3D(
                    abm::util::toCartesianCoordinates(*it2) - abm::util::toCartesianCoordinates(*it));
            double distance = dConnect.getMagnitude();
            if (distance < radiusPoresOfKohn + radiusAlvEpithTypeOne && distance > 0.001) {
                currentEpithelialNeigbourCellNode = xmlTags->addChildToNode(currentEpithelialCellNode, "PoKNC");
                Coordinate3D neighbourPore = abm::util::toCartesianCoordinates(*it2);
                std::ostringstream ssRadiusPoK;
                ssRadiusPoK << radiusPoresOfKohn;
                xmlTags->addDataFieldToNode(currentEpithelialNeigbourCellNode, "radius", "discrete", "double",
                                            ssRadiusPoK.str());
                xmlTags->addDataFieldToNode(currentEpithelialNeigbourCellNode, neighbourPore);
            }

        }
    }

    //include type II cells into the output
    for (auto it = alvEpithTypeTwo.begin(); it != alvEpithTypeTwo.end(); it++) {

        currentEpithelialCellNode = xmlTags->addChildToNode(sphereSiteNode, "AET2");
        Coordinate3D curPos = abm::util::toCartesianCoordinates(*it);

        std::ostringstream ssLengthEpT;
        ssLengthEpT << lengthAlvEpithTypeTwo;

        xmlTags->addDataFieldToNode(currentEpithelialCellNode, "length", "discrete", "double", ssLengthEpT.str());
        xmlTags->addDataFieldToNode(currentEpithelialCellNode, curPos);
    }

    //include Pores of Kohn into the output
    for (auto it = poresOfKohn.begin(); it != poresOfKohn.end(); it++) {

        currentEpithelialCellNode = xmlTags->addChildToNode(sphereSiteNode, "PoK");
        Coordinate3D curPos = abm::util::toCartesianCoordinates(*it);

        std::ostringstream ssRadiusPoK;
        ssRadiusPoK << radiusPoresOfKohn;

        xmlTags->addDataFieldToNode(currentEpithelialCellNode, "radius", "discrete", "double", ssRadiusPoK.str());
        xmlTags->addDataFieldToNode(currentEpithelialCellNode, curPos);
    }

}

void AlveoleSite::includeAlveolarEpitheliumType1() {
    radiusAlvEpithTypeOne = r0AEC1;
    int countAEC1 = 0;
    double phi, theta;
    int objectsPerRow = opR;
    SphericCoordinate3D sC3d;
    for (theta = thetaLowerBound; theta <= M_PI; theta += M_PI / (0.5 * objectsPerRow)) {
        double numberOfObjectsPerPhiRow = std::max(floor(1.0 * objectsPerRow * sin(theta)), 1.0);
        for (phi = 0; phi < 2 * M_PI; phi += 2 * M_PI / numberOfObjectsPerPhiRow) {
            sC3d = SphericCoordinate3D{radius + thicknessOfBorder / 2.0, theta, phi};
            alvEpithTypeOne.push_back(sC3d);
            countAEC1++;
        }
    }
    std::vector<SphericCoordinate3D>::iterator it;
    unsigned int m = 0, n = 0;
    for (it = alvEpithTypeOne.begin(); it != alvEpithTypeOne.end(); it++) {
        auto it2 = alvEpithTypeOne.begin();
        advance(it2, n + 1);
        m = n + 1;
        for (; it2 != alvEpithTypeOne.end(); it2++) {
            Coordinate3D dConnect = Coordinate3D(
                    abm::util::toCartesianCoordinates(*it2) - abm::util::toCartesianCoordinates(*it));
            double distance = dConnect.getMagnitude();
            if (distance < 1.5 * radiusAlvEpithTypeOne && distance > 0.001) {
                std::pair<int, int> curPair(m, n);
                pairCells.push_back(curPair);
            }
            m++;
        }
        n++;
    }
}

void AlveoleSite::includePoresOfKohnHuman() {
    int m, n, iterations = 0;
    // Place the Pores of Kohn at the intersections between neighbouring AECT1
    for (int j = 0; j < noOfPoK; j++) {
        SphericCoordinate3D sphericPositionPoK;
        std::vector<std::pair<int, int> >::iterator itPairs;
        iterations = 0;
        double mindist = 1000000;
        bool condition1Theta, condition2Theta, MaxItNotReached, MindistReached;
        do {
            if (!pairCells.empty()) {
                unsigned int rand = random_generator_->generateInt(pairCells.size() - 1);
                itPairs = pairCells.begin();
                advance(itPairs, rand);
                m = itPairs->first;
                n = itPairs->second;
                Coordinate3D c3d{abm::util::toCartesianCoordinates(alvEpithTypeOne[n]) -
                                 abm::util::toCartesianCoordinates(alvEpithTypeOne[m])};
                c3d *= 0.5;
                Coordinate3D positionPoK = abm::util::toCartesianCoordinates(alvEpithTypeOne[m]);
                positionPoK += c3d;
                sphericPositionPoK = abm::util::toSphericCoordinates(positionPoK);
                sphericPositionPoK.r = radius;
                mindist = minDistanceToPoK(sphericPositionPoK);
                iterations++;
            } else {
                break;
            }
            condition1Theta = sphericPositionPoK.theta < M_PI - (30.0 / 116.5);
            condition2Theta = sphericPositionPoK.theta > thetaLowerBound + 0.35;
            MaxItNotReached = iterations < 1000;
            MindistReached = mindist > radiusPoresOfKohn * 3;
        } while (!(condition1Theta && condition2Theta && MaxItNotReached && MindistReached));
        poresOfKohn.push_back(sphericPositionPoK);
        if (!pairCells.empty()) {
            pairCells.erase(itPairs);
        }
    }
}

void AlveoleSite::includeAlveolarEpitheliumType2Human() {
    unsigned int m = 0, n = 0, iterations;
    // Place the AEC2 at the intersections between neighbouring AEC1
    for (int j = 0; j < noOfAEC2; j++) {
        std::vector<std::pair<int, int> >::iterator itPairs;
        SphericCoordinate3D sphericPositionAECT2{};
        iterations = 0;
        do {
            if (!pairCells.empty()) {
                unsigned int rand = random_generator_->generateInt(pairCells.size() - 1);
                itPairs = pairCells.begin();
                advance(itPairs, rand);
                m = itPairs->first;
                n = itPairs->second;
                Coordinate3D c3d{abm::util::toCartesianCoordinates(alvEpithTypeOne[n]) -
                                 abm::util::toCartesianCoordinates(alvEpithTypeOne[m])};
                c3d *= 0.5;
                Coordinate3D positionAECT2 = abm::util::toCartesianCoordinates(alvEpithTypeOne[m]);
                positionAECT2 += c3d;
                sphericPositionAECT2 = abm::util::toSphericCoordinates(positionAECT2);
                sphericPositionAECT2.r = radius + thicknessOfBorder / 2.0;
                iterations++;
            } else {
                break;
            }
        } while (iterations < 1000);

        if (!pairCells.empty()) {
            pairCells.erase(itPairs);
        } else {
            break;
        }

        if (minDistanceToPoK(sphericPositionAECT2) < radiusPoresOfKohn * 2) {
            j--;
        } else {
            alvEpithTypeTwo.push_back(sphericPositionAECT2);
        }

    }
}

void AlveoleSite::calculateCrossPoints() {
    crossPoints.clear();
    while (crossPoints.size() < 50000) {
        Coordinate3D pos = getRandomMinDistanceToBoundaryPosition(5);
        std::vector<double> distances;
        for (size_t j = 0; j < alvEpithTypeOne.size(); j++) {
            distances.push_back(
                    alvEpithTypeOne.at(j).calculateSphericalDistance(abm::util::toSphericCoordinates(pos)));
        }
        std::sort(distances.begin(), distances.end(), [](double u, double v) -> bool { return u > v; });
        for (size_t j = 0; j < distances.size(); j++) {
        }
        if (abs(distances.at(distances.size() - 1) - distances.at(distances.size() - 2)) < 3) {
            crossPoints.push_back(pos);
        }
    }
}

void AlveoleSite::includePoresOfKohnMouse() {
    calculateCrossPoints();
    double mindist{};
    Coordinate3D posPOK{};
    int iterations = 0;
    // Place the Pores of Kohn at the intersections between neighbouring AECT1
    for (int j = 0; j < noOfPoK; j++) {
        do {
            mindist = 1000000;
            // Potential PoK positions were determined in crosspoints, now we pick
            int a = random_generator_->generateInt(crossPoints.size() - 1);
            posPOK = crossPoints.at(a);
            for (const auto &pore : poresOfKohn) {
                double dist = posPOK.calculateEuclidianDistance(abm::util::toCartesianCoordinates(pore));
                if (dist < mindist) mindist = dist;
            }
            iterations++;
        } while (mindist < 6 && iterations < 1000); //fpt measurement condition
        if (iterations < 1000) {
            poresOfKohn.push_back(abm::util::toSphericCoordinates(posPOK));
        } else {
            ERROR_STDERR("could not find a position for the pore of Kohn!!!");
        }
    }
}

void AlveoleSite::includeAlveolarEpitheliumType2Mouse() {
    double mindistPOK{}, mindistAEC2{};
    Coordinate3D posAEC2;
    int iterations = 0;
    // Place the AEC2 at the intersections between neighbouring AEC1
    for (int j = 0; j < noOfAEC2; j++) {
        do {
            mindistPOK = 1000000;
            mindistAEC2 = 1000000;
            // Potential AEC2 positions were determined in crosspoints, now we pick
            int a = random_generator_->generateInt(crossPoints.size() - 1);
            posAEC2 = crossPoints.at(a);
            for (size_t k = 0; k < poresOfKohn.size(); k++) {
                double dist = posAEC2.calculateEuclidianDistance(abm::util::toCartesianCoordinates(poresOfKohn.at(k)));
                if (dist < mindistPOK) mindistPOK = dist;
            }
            for (const auto &aec2 : alvEpithTypeTwo) {
                double dist = posAEC2.calculateEuclidianDistance(abm::util::toCartesianCoordinates(aec2));
                if (dist < mindistAEC2) mindistAEC2 = dist;
            }
            iterations++;
        } while ((mindistPOK < 10 || mindistAEC2 < 15) && iterations < 1000); //fpt measurement condition
        if (iterations < 1000) {
            alvEpithTypeTwo.push_back(abm::util::toSphericCoordinates(posAEC2));
        } else {
            ERROR_STDERR("could not find a position for the AEC2!!!");
        }
    }
}

Coordinate3D AlveoleSite::generateRandomDirectionVector(Coordinate3D position, double length) {
    double x, y, z, subst, r, phi, theta, u, dtheta, alpha, vix, viy, viz, radiusOrbit;
    SphericCoordinate3D currentPos = abm::util::toSphericCoordinates(position - centerOfSite);

    switch (dimensions) {
        case 2:
            radiusOrbit = currentPos.r;
            dtheta = length / radiusOrbit; //in rad
            alpha = random_generator_->generateDouble(M_PI * 2.0); //direction of the vector
            alpha2dTurningAngle = alpha;

            vix = radiusOrbit * sin(dtheta) * cos(alpha);
            viy = radiusOrbit * sin(dtheta) * sin(alpha);
            viz = radiusOrbit * (cos(dtheta) - 1);
            phi = currentPos.phi;
            theta = currentPos.theta;

            x = cos(phi) * (cos(theta) * vix + sin(theta) * viz) - sin(phi) * viy;
            y = sin(phi) * (cos(theta) * vix + sin(theta) * viz) + cos(phi) * viy;
            z = cos(theta) * viz - sin(theta) * vix;

            break;

        default:
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

Coordinate3D AlveoleSite::generateBiasedRandomDirectionVector(Agent *agent, Coordinate3D position, double length) {

    Coordinate3D result{}, randomPart{}, directedPart{}, intermediatePositions{};
    intermediatePositions = position;

    double p = agent->getFeatureValueByName("migration-bias-probability"); // probability to move up the gradient
    double alpha = agent->getFeatureValueByName("gradient-direction"); // direction of the gradient as sensed by the

    // Biased Persistent Random Walk with mixed parts of random and directed walk
    if (random_generator_->generateInt(1) > 0) {
        randomPart = generateRandomDirectionVector(position, (1.0 - p) * length);
        intermediatePositions += randomPart;
        directedPart = generateDirectedVector(intermediatePositions, alpha, p * length);
        intermediatePositions += directedPart;
    } else {
        directedPart = generateDirectedVector(intermediatePositions, alpha, p * length);;
        intermediatePositions += directedPart;
        randomPart = generateRandomDirectionVector(intermediatePositions, (1.0 - p) * length);
        intermediatePositions += randomPart;
    }

    result = generateDirectedVector(position, abm::util::toSphericCoordinates(intermediatePositions - centerOfSite),
                                    length);
    return result;
}

Coordinate3D AlveoleSite::generatePersistentDirectionVector(Coordinate3D position,
                                                            double length,
                                                            Coordinate3D prevVector,
                                                            double previousAlpha) {
    double x, y, z, subst, r, phi, theta, u, dtheta, alpha, vix, viy, viz, radiusOrbit;
    SphericCoordinate3D currentPos = abm::util::toSphericCoordinates(position - centerOfSite);

    switch (dimensions) {
        case 2:
            radiusOrbit = currentPos.r;
            dtheta = length / radiusOrbit; //in rad
            alpha = previousAlpha; //direction of the vector

            vix = radiusOrbit * sin(dtheta) * cos(alpha);
            viy = radiusOrbit * sin(dtheta) * sin(alpha);
            viz = radiusOrbit * (cos(dtheta) - 1);
            phi = currentPos.phi;
            theta = currentPos.theta;

            x = cos(phi) * (cos(theta) * vix + sin(theta) * viz) - sin(phi) * viy;
            y = sin(phi) * (cos(theta) * vix + sin(theta) * viz) + cos(phi) * viy;
            z = cos(theta) * viz - sin(theta) * vix;

            break;

        default:
            x = prevVector.x;
            y = prevVector.y;
            z = prevVector.z;

            break;
    }

    return Coordinate3D{x, y, z};
}

Coordinate3D AlveoleSite::generateDirectedVector(Coordinate3D position, SphericCoordinate3D posOfGoal, double length) {
    double x = 0, y = 0, z = 0, phi = 0, theta = 0, dtheta = 0, alpha = 0, vix = 0, viy = 0, viz = 0, radiusOrbit = 0;
    SphericCoordinate3D currentPos = abm::util::toSphericCoordinates(position - centerOfSite);

    switch (dimensions) {
        case 2:
            radiusOrbit = currentPos.r;
            dtheta = length / radiusOrbit; //in rad
            alpha = retrieveDirectionAngleAlpha(currentPos, posOfGoal); //direction of the vector
            alpha2dTurningAngle = alpha;

            vix = radiusOrbit * sin(dtheta) * cos(alpha);
            viy = radiusOrbit * sin(dtheta) * sin(alpha);
            viz = radiusOrbit * (cos(dtheta) - 1);
            phi = currentPos.phi;
            theta = currentPos.theta;

            x = cos(phi) * (cos(theta) * vix + sin(theta) * viz) - sin(phi) * viy;
            y = sin(phi) * (cos(theta) * vix + sin(theta) * viz) + cos(phi) * viy;
            z = cos(theta) * viz - sin(theta) * vix;

            break;

        default:
            x = 0;
            y = 0;
            z = 0;

            break;
    }

    return Coordinate3D{x, y, z};
}

Coordinate3D AlveoleSite::generateDirectedVector(Coordinate3D position, double alpha, double length) {
    double x = 0, y = 0, z = 0, phi = 0, theta = 0, dtheta = 0, vix = 0, viy = 0, viz = 0, radiusOrbit = 0;
    SphericCoordinate3D currentPos = abm::util::toSphericCoordinates(position - centerOfSite);

    switch (dimensions) {
        case 2:
            radiusOrbit = currentPos.r;
            dtheta = length / radiusOrbit; //in rad
            alpha2dTurningAngle = alpha;

            vix = radiusOrbit * sin(dtheta) * cos(alpha);
            viy = radiusOrbit * sin(dtheta) * sin(alpha);
            viz = radiusOrbit * (cos(dtheta) - 1);
            phi = currentPos.phi;
            theta = currentPos.theta;

            x = cos(phi) * (cos(theta) * vix + sin(theta) * viz) - sin(phi) * viy;
            y = sin(phi) * (cos(theta) * vix + sin(theta) * viz) + cos(phi) * viy;
            z = cos(theta) * viz - sin(theta) * vix;

            break;

        default:
            x = 0;
            y = 0;
            z = 0;

            break;
    }

    return Coordinate3D{x, y, z};
}

bool AlveoleSite::overAECT1(SphericCoordinate3D posConidia) {
    posObstacle = posConidia;
    double minDistance = 1000000;
    std::vector<SphericCoordinate3D>::iterator it2;
    for (it2 = alvEpithTypeOne.begin(); it2 != alvEpithTypeOne.end(); it2++) {
        if (it2->calculateSphericalDistance(posObstacle) < minDistance) {
            minDistance = it2->calculateSphericalDistance(posObstacle);
            cellOfObstacle = *it2;
        }
    }
    obstacleIsOnType1 = true;

//    check if obstacle is on type II
    for (it2 = alvEpithTypeTwo.begin(); it2 != alvEpithTypeTwo.end(); it2++) {
        double arcOfIntervalTheta = 0.5 * lengthAlvEpithTypeTwo / (radius);
        double thetaObstacle = posObstacle.theta;
        double thetaAEC2 = it2->theta;
        //first condition: the spore has to be in a defined theta range [\theta-d\theta, \theta+d\theta]
        if ((thetaObstacle <= thetaAEC2 + arcOfIntervalTheta) && (thetaObstacle >= thetaAEC2 - arcOfIntervalTheta)) {
            double d = it2->calculateEuclidianDistance(posObstacle);
            double dmax =
                    sqrt(lengthAlvEpithTypeTwo * lengthAlvEpithTypeTwo / 4.0 +
                         (thetaObstacle - thetaAEC2) * (thetaObstacle - thetaAEC2) *
                         radius * radius);
            // second condition: the spore should have a maximum distance from
            if (d <= dmax) {
                cellOfObstacle = *it2;
                obstacleIsOnType1 = false;
                break;
            }
        }
    }
    return obstacleIsOnType1;
}

double AlveoleSite::retrieveDirectionAngleAlpha(SphericCoordinate3D ownPos, SphericCoordinate3D goalPos) {
    double alpha = 0;
    double dphi = ownPos.phi - goalPos.phi;
    double dphiAbs = fabs(dphi);
    double c = ownPos.theta;
    double b = goalPos.theta;
    double acosarg = sin(ownPos.theta) * sin(goalPos.theta) * cos(dphi) + cos(ownPos.theta) * cos(goalPos.theta);
    if (abs(acosarg) > 1) INFO_STDOUT("WARNING: ACOS argument might be greater than zero");
    double a = (abs(acosarg) < 1) ? acos(acosarg) : 0;
    if (c == 0) {
        alpha = goalPos.phi;
    } else {
        if (b == 0) {
            alpha = M_PI;
        } else {
            //compute beta by cases
            double beta;
            if (ownPos.phi == goalPos.phi) {
                if (ownPos.theta > goalPos.theta) {
                    beta = 0;
                } else {
                    beta = M_PI;
                }
            } else {
                acosarg = (cos(b) - cos(a) * cos(c)) / (sin(a) * sin(c));
                if (abs(acosarg) > 1) INFO_STDOUT("WARNING: ACOS argument might be greater than zero");
                beta = (abs(acosarg) < 1) ? acos(acosarg) : 0;
            }
            //compute alpha by cases
            if ((ownPos.phi <= goalPos.phi && dphiAbs <= M_PI) || (ownPos.phi > goalPos.phi && dphiAbs > M_PI)) {
                alpha = M_PI - beta;
            } else {
                alpha = M_PI + beta;
            }
        }
    }
    return alpha;
}

bool AlveoleSite::onAECTObstacleCell(Coordinate3D position) {
    bool onAECTObstacle = false;
    if (obstacleIsOnType1) {
        double dPosGoal = abm::util::toSphericCoordinates(position).calculateSphericalDistance(cellOfObstacle);
        if (dPosGoal < 2 * radiusAlvEpithTypeOne) {
            // Check for even closer AEC1
            std::vector<SphericCoordinate3D>::iterator it2;
            for (it2 = alvEpithTypeOne.begin(); it2 != alvEpithTypeOne.end(); it2++) {
                if (abm::util::toSphericCoordinates(position).calculateSphericalDistance((*it2)) < dPosGoal) {
                    break;
                }
            }
            if (it2 != alvEpithTypeOne.end()) {
                onAECTObstacle = false;
            } else {
                onAECTObstacle = true;
            }
        } else {
            onAECTObstacle = false;
        }
    } else {
        double arcOfIntervalTheta = 0.5 * lengthAlvEpithTypeTwo / (radius);
        double arcOfIntervalPhi = 0.5 * lengthAlvEpithTypeTwo / (radius * sin(cellOfObstacle.theta));
        double phiPos = abm::util::toSphericCoordinates(position).phi;
        double thetaPos = abm::util::toSphericCoordinates(position).theta;
        if (phiPos < cellOfObstacle.phi + arcOfIntervalPhi && phiPos > cellOfObstacle.phi - arcOfIntervalPhi &&
            thetaPos < cellOfObstacle.theta + arcOfIntervalTheta &&
            thetaPos > cellOfObstacle.theta - arcOfIntervalTheta) {
            onAECTObstacle = true;
        }
    }
    return onAECTObstacle;
}

Coordinate3D AlveoleSite::getRandomPosition() {
    Coordinate3D randPos{};
    do {
        randPos = this->SphereSite::getRandomPosition();
    } while (!containsPosition(randPos));
    return randPos;
}

Coordinate3D AlveoleSite::getRandomBoundaryPoint() {
    Coordinate3D boundaryPoint{};
    // Decide whether to use one of pores of Kohn or the alveolar entrance ring
    double lengthOfPoKLineElements = 2 * M_PI * radiusPoresOfKohn * noOfPoK;
    double lengthOfAERLineElements = 2 * M_PI * radius * sin(thetaLowerBound);
    do {
        double decisionRand = random_generator_->generateDouble(lengthOfPoKLineElements + lengthOfAERLineElements);
        if (decisionRand < lengthOfPoKLineElements) {
            // Use a pore of Kohn as boundary point
            unsigned int entrancePoKindex = random_generator_->generateInt(noOfPoK - 1);
            boundaryPoint = abm::util::toCartesianCoordinates(poresOfKohn[entrancePoKindex]);
            Coordinate3D shiftFromPoKCenter = generateRandomDirectionVector(
                    boundaryPoint, radiusPoresOfKohn * thetaLowerBound);
            boundaryPoint += shiftFromPoKCenter;
            DEBUG_STDOUT("entering an AM at PoK");
        } else {
            // Use the alveolar entrance ring as boundary point
            double r, phi, theta;
            SphericCoordinate3D sc3d{};
            boundaryPoint = getRandomPosition();
            sc3d = abm::util::toSphericCoordinates(boundaryPoint - centerOfSite);
            r = sc3d.r;
            phi = sc3d.phi;
            theta = thetaLowerBound + 0.001;
            sc3d = SphericCoordinate3D{r, theta, phi};
            boundaryPoint = abm::util::toCartesianCoordinates(sc3d);
        }
    } while (!containsPosition(boundaryPoint));

    // Generate a boundary vector in order to keep the agent within the system at least for the next step
    Coordinate3D nextPos = Coordinate3D();
    do {
        boundary_input_vector_ = generateRandomDirectionVector(boundaryPoint, 1);
        nextPos = boundaryPoint;
        nextPos += boundary_input_vector_;
    } while (!containsPosition(nextPos));

    return boundaryPoint;
}

double AlveoleSite::getDistanceFromBoundary(Coordinate3D position) {
    double minDistance = 1000000;

    if (containsPosition(position)) {
        //check alveolar entrance ring
        SphericCoordinate3D posSpheric = abm::util::toSphericCoordinates(position);
        double dtheta = posSpheric.theta - thetaLowerBound;
        minDistance = dtheta * posSpheric.r;

        // Check for pores of Kohn
        std::vector<SphericCoordinate3D>::iterator it;
        double distance;
        for (it = poresOfKohn.begin(); it != poresOfKohn.end(); it++) {
            SphericCoordinate3D poreOfKohnPos = *it;
            distance = posSpheric.calculateSphericalDistance(poreOfKohnPos);
            if (distance - radiusPoresOfKohn < minDistance) {
                minDistance = distance - radiusPoresOfKohn;
            }
        }
    }
    if (minDistance < 0) minDistance = 0;

    return minDistance;
}

Coordinate3D AlveoleSite::getRandomMinDistanceToBoundaryPosition(double minDistanceToBoundary) {
    Coordinate3D position{};
    do {
        position = getRandomPosition();
    } while (getDistanceFromBoundary(position) < minDistanceToBoundary);

    return position;
}

double AlveoleSite::minDistanceToPoK(const SphericCoordinate3D &sc3d) {
    double minDistance = 1000000, distance; //µm

    std::vector<SphericCoordinate3D>::iterator it2;
    for (it2 = poresOfKohn.begin(); it2 != poresOfKohn.end(); it2++) {
        distance = it2->calculateSphericalDistance(sc3d);
        if (distance < minDistance) {
            minDistance = distance;
        }
    }
    return minDistance;
}

double AlveoleSite::getFeatureValueByName(std::string name) {
    double value = 0;
    if (name == "surfactantThickness") {
        value = surfactantThickness;
    }
    if (name == "radius") {
        value = radius;
    }
    if (name == "radius") {
        value = radius;
    }
    if (name == "radiusPoresOfKohn") {
        value = radiusPoresOfKohn;
    }
    if (name == "radiusAlvEpithTypeOne") {
        value = radiusAlvEpithTypeOne;
    }
    if (name == "lengthAlvEpithTypeTwo") {
        value = lengthAlvEpithTypeTwo;
    }
    return value;
}

void AlveoleSite::initializeAgents(const abm::util::SimulationParameters::AgentManagerParameters &parameters,
                                   const std::string &input_dir,
                                   double current_time,
                                   double time_delta) {

    agent_manager_->setLastConidiaChange(0.0);
    std::vector<Coordinate3D> AMpos;

    // Loop over agent types
    for (const auto &agent_parameters: parameters.agents) {
        auto initDistribution = agent_parameters->initial_distribution;
        std::string inputDistributionPath = boost::filesystem::path(input_dir).append(
                agent_parameters->input_distribution_path).string();
        const auto name = agent_parameters->type;
        int numberOfAgents = static_cast<int>(agent_parameters->number);

        if (name == "Macrophage") {
            std::stringstream entirePath;
            if (numberOfAgents < 0 and organism == 1) {
                entirePath << inputDistributionPath << "AM4-375" << "/";
            } else {
                int pathAM = (organism == 1) ? numberOfAgents : static_cast<int>(agent_parameters->number * 10);
                entirePath << inputDistributionPath << "AM" << std::to_string(pathAM) << "/";
            }
            // baseInputRates are the mean value over all calibrated lambda values / AM numbers for each system and result from previous analysis
            double baseInputRate = (organism == 1) ? 0.011863324353554 : 0.051869337535733;
            if (initDistribution == 1 && !abm::util::folderExists(entirePath.str())) {
                initDistribution = 0;
                // actual input rate correlates linearly with AM number
                inputRate = baseInputRate * agent_parameters->number;
                ERROR_STDERR(
                        "WARNING: For given AM number, input distribution " + entirePath.str() + " does not exist. \nChoose randomly and adjust lambda input rate.");
            } else if (initDistribution == 0 && inputRate == 0.0) {
                // actual input rate correlates linearly with AM number
                inputRate = baseInputRate * agent_parameters->number;
                ERROR_STDERR("WARNING: Input rate is 0.0. Change it to constant AM inputs.");
            }

            // Determine the numberOfAgents to put into the system
            if (initDistribution == 0 && agent_parameters->binomial_distribution.activated) {
                int curNo{};
                double p = agent_parameters->binomial_distribution.p; // p ~ 1/Total number of Alveoli in lung
                auto n = agent_parameters->binomial_distribution.n; // n ~ Total number of AM in lung
                double pRand{};
                do {
                    curNo = static_cast<int>(getRandomGenerator()->generateInt(30));
                    pRand = getRandomGenerator()->generateDouble();
                } while (pRand > Algorithms::bernoulliProbability(n, curNo, p));
                numberOfAgents = curNo;
            }

            // Load AM distribution from input files and corresponding lambda value
            if (initDistribution == 1) {
                DEBUG_STDOUT("Choose AM positions from distributions.");
                inputDistributionPath = entirePath.str();
                // Choose input rate from calibrated lambda for the specific AM number
                std::stringstream pathToLambdaValue;
                pathToLambdaValue << inputDistributionPath << "lambdaInputRate.txt";
                inputRate = abm::util::readLambdaValueFromFile(pathToLambdaValue.str());
                std::vector<std::string> fileNames = abm::util::getFileNamesFromDirectory(inputDistributionPath, "POK");
                bool AMinside = false;
                while (!AMinside) {
                    AMpos.clear();
                    // Pick AM positions randomly from input file
                    int pick = getRandomGenerator()->generateInt(0, fileNames.size() - 1);
                    std::stringstream pathToFile;
                    pathToFile << inputDistributionPath << fileNames[pick];
                    abm::util::read3DCoordinatesFromFile(AMpos, pathToFile.str());
                    DEBUG_STDOUT("Insert " + std::to_string(AMpos.size()) + " AMs from " + pathToFile.str());
                    // Check if all AM are inside the alveolus, it may happen that they are placed over a PoK
                    AMinside = true;
                    for (auto &AMpo : AMpos) {
                        if (!containsPosition(AMpo)) {
                            AMinside = false;
                        }
                    }
                }
                numberOfAgents = AMpos.size();
            }
        }

        // Create agents
        for (auto agent_pos = 0; agent_pos < numberOfAgents; agent_pos++) {
            Coordinate3D initialPosition{};
            std::shared_ptr<Agent> agent{};

            switch (initDistribution) {
                case 0:
                    initialPosition = getRandomPosition();
                    break;
                case 1:
                    // AM is placed according to a distribution which is loaded from file
                    if (name == "Macrophage") {
                        DEBUG_STDOUT("Insert AM" << agent_pos << " at position " << AMpos.at(agent_pos).printCoordinates()
                                            << ", " << containsPosition(AMpos.at(agent_pos)) << ", "
                                            << containsPosition(AMpos.at(agent_pos)));
                        initialPosition = AMpos.at(agent_pos);
                    } else {
                        initialPosition = getRandomPosition();
                    }
                    break;
                default:
                    initialPosition = getRandomPosition();
                    break;
            }
            agent = agent_manager_->emplace_back(CellFactory::createCell(name, std::make_unique<Coordinate3D>(initialPosition),
                                                                         agent_manager_->generateNewID(),this,
                                                                         time_delta, current_time));

            if (name == "AspergillusFumigatus") agent_manager_->addConidiaToList(agent.get());
        }

    }
}
