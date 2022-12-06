//  Copyright by Christoph Saffer, Paul Rudolph, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#include <algorithm>
#include <map>

#include "simulation/ParticleManager.h"
#include "simulation/Site.h"
#include "analyser/InSituMeasurements.h"
#include "utils/macros.h"
#include "simulation/AgentManager.h"


ParticleManager::ParticleManager(Site *site) {
    this->site = site;
    noOfParticles = 0;
    idHandling = 0;
    allowHigherDT = false;
    conidiumParticle = 0;
}

void ParticleManager::initializeParticles(Site *site,
                                          const abm::util::SimulationParameters::ParticleManagerParameters &parameters,
                                          const std::string &input_dir) {

    //initDistribution: 0-random, 1-every agent in one place at beginning
    dc = parameters.diffusion_constant;
    particleInputDelauneyFile = parameters.particle_delauney_input_file;
    drawIsolines = parameters.draw_isolines;

    particleSecretionMoleculePerCellMin = parameters.molecule_secretion_per_cell; //6000.0;

    particleBalloonList = std::make_unique<StaticBalloonList>(10.61, site->getLowerLimits(), site->getUpperLimits());
    particleBalloonList->setThreshold(10.6);
    triangulationFromDirectInput(site, parameters, input_dir);

}

void ParticleManager::triangulationFromDirectInput(Site *site,
                                                   const abm::util::SimulationParameters::ParticleManagerParameters &parameters,
                                                   const std::string &input_dir) {

    auto dc = parameters.diffusion_constant;
    double R = site->getRadius();

    std::vector<Coordinate3D> coordinates;
    InputConfiguration icParticles = InputConfiguration(boost::filesystem::path(input_dir).append(particleInputDelauneyFile).string());
    XMLNode particlesForInput = icParticles.getRootNode().getChildNode("Particles");

    // extract particles first
    for (int i = 0; i < particlesForInput.nChildNode(); i++) {
        XMLNode curParticle = particlesForInput.getChildNode(i);
        std::string name = curParticle.getName();
        if (name.compare("Particle") == 0) {

            Coordinate3D position = InputConfiguration::getCoordinateDataFieldValueByName(&curParticle, "position");
            double area = InputConfiguration::getDoubleDataFieldValueByName(&curParticle, "area");
            double conc = InputConfiguration::getDoubleDataFieldValueByName(&curParticle, "concentration");
            coordinates.push_back(position);

            auto particle = createParticle(site, position, conc, dc);
            particle->setArea(area);
        }
    }

    // extract neighbours of each particle and set that shit!
    for (int i = 0; i < particlesForInput.nChildNode(); i++) {
        XMLNode curParticle = particlesForInput.getChildNode(i);
        std::string name = curParticle.getName();
        if (name.compare("Particle") == 0) {
            unsigned int id = InputConfiguration::getUIntDataFieldValueByName(&curParticle, "id");
            Particle *p = particleById[id];
            XMLNode interactionPartners = curParticle.getChildNode("Interactions");
            for (int j = 0; j < interactionPartners.nChildNode(); j++) {
                XMLNode curInteractionPartner = interactionPartners.getChildNode(j);
                std::string partnerNodeName = curInteractionPartner.getName();
                if (partnerNodeName == "Particle") {
                    unsigned int partnerId = InputConfiguration::getUIntDataFieldValueByName(&curInteractionPartner,"id");
                    double contactArea = InputConfiguration::getDoubleDataFieldValueByName(&curInteractionPartner, "contact");
                    Particle *pNeighbour = particleById[partnerId];
                    p->addNeighbour(pNeighbour, contactArea);
                }
            }
        }
    }
    computeMinTimestepDistribution();

    // draw Isolines in Visualisation
    if (drawIsolines) extractTriangles();

    // set particles as inside or outside the environment
    int inside = 0, boundary = 0, outside = 0;

    auto itP = allParticles.begin();
    while (itP != allParticles.end()) {

        (*itP)->setLabel();

        if ((*itP)->getIsInSite()) {
            inside++;
        } else {
            if ((*itP)->getIsAtBoundary()) {
                boundary++;
            } else {
                outside++;
            }
        }

        itP++;
    }

    DEBUG_STDOUT("particle statistics: inside:" + std::to_string(inside) +
                 " boundary:" + std::to_string(boundary) +
                 " outside:" + std::to_string(outside) +
                 " sum:" + std::to_string(inside + outside + boundary));

    std::shared_ptr<Particle> pAtMPI = nullptr;
    double areaFull = 4 * M_PI * R * R;
    double areaMeasure = 0;

    std::vector<double> contactAreas, distances;
    while (itP != allParticles.end()) {
        areaMeasure += (*itP)->getArea();
        contactAreas = (*itP)->getParticleNeighbourList()->getContactAreas();
        distances = (*itP)->getParticleNeighbourList()->getDistances();
        auto itPN = (*itP)->getParticleNeighbourList()->getNeighbours().begin();
        int idx = 0;
        while (itPN != (*itP)->getParticleNeighbourList()->getNeighbours().end()) {
            if (distances[idx] == 0) {
                DEBUG_STDOUT(
                        "distance " + std::to_string(distances[idx]) +
                        " for ids:" + std::to_string((*itP)->getId()) + " " +
                        std::to_string((*itPN)->getId()) + " " +
                        std::to_string((*itP)->getPosition().calculateEuclidianDistance(
                                (*itPN)->getPosition())) + " ");
            }
            itPN++;
            idx++;
        }
        itP++;
    }
    DEBUG_STDOUT("area correct:" + std::to_string(areaFull) +
                 " measured/computed from cells:" + std::to_string(areaMeasure) +
                 " rel.Error:" + std::to_string(1 - areaMeasure / areaFull));
}

std::shared_ptr<Particle> ParticleManager::getParticleByPosition(Coordinate3D pos) {
    std::shared_ptr<Particle> p = nullptr;
    for (auto itP: allParticles){
        double dist = itP->getPosition().calculateEuclidianDistance(pos);
        if (dist < 0.001) {
            p = itP;
            break;
        }
    }
    return p;
}

void ParticleManager::computeMinTimestepDistribution() {
    double minTimestep = 1000000, timestep{};
    for (auto itP: allParticles){
        timestep = itP->estimateLowestTimestep();
        if (timestep < minTimestep)
            minTimestep = timestep;
    }
    INFO_STDOUT("maximum timestep that is possible for stability reasons: " + std::to_string(minTimestep));
}

void ParticleManager::inputOfParticles(double time_delta) {

    if (clean_chemotaxis && !site->getLargeTimestepActive()) {
        DEBUG_STDOUT("Start chemotaxis cleanup");
        cleanUpAllParticles();
        clean_chemotaxis = false;
    }
    insertConcentrationAtArea(site, time_delta);
}

void ParticleManager::insertConcentrationAtArea(Site *site, double time_delta) {
    if (aecParticles.empty() && particleSecretionMoleculePerCellMin > 0) {
        sumAreaAECParticles = 0;
        aecParticlesCells.clear();
        aecParticlesCells.clear();
        for (int i = 0; i < 100; i++) {
            sumAreaAEcParticlesCells[i] = 0;
            aecSecretionratePerGrid[i] = 0;
        }
        std::vector<Agent *> allConidia = site->getAgentManager()->getAllConidia();
        for (size_t i = 0; i < allConidia.size(); i++) {
            bool overAECT1;
            if (!allConidia.at(i)->isDeleted()) {
                overAECT1 = site->overAECT1(abm::util::toSphericCoordinates(allConidia.at(i)->getPosition()));
                SphericCoordinate3D posObstacleAEC = abm::util::toSphericCoordinates(allConidia.at(i)->getPosition());
                std::vector<unsigned int> potentialAECParticles;
                int ID;
                if (overAECT1) {
                    ID = get_closest_AEC_ID(allConidia.at(i)->getPosition(), 1);
                    INFO_STDOUT("Found Conidia " + std::to_string(i) + " at position " << allConidia.at(i)->getPosition().printCoordinates() << " over AECI " + std::to_string(ID));
                    particleBalloonList->setThreshold(48.0);
                } else {
                    ID = get_closest_AEC_ID(allConidia.at(i)->getPosition(), 2);
                    INFO_STDOUT("Found Conidia " + std::to_string(i) + " at position " << allConidia.at(i)->getPosition().printCoordinates() << " over AECII " + std::to_string(ID));
                    particleBalloonList->setThreshold(6.61);
                }
                particleBalloonList->getInteractions(abm::util::toCartesianCoordinates(posObstacleAEC), potentialAECParticles);

                auto itP = potentialAECParticles.begin();
                while (itP != potentialAECParticles.end()) {
                    auto currentParticle = allParticles[(*itP)];
                    if (currentParticle->getIsInSite()) {
                        if (site->onAECTObstacleCell(currentParticle->getPosition())){
                            if (std::find(aecParticles.begin(), aecParticles.end(), currentParticle) ==
                                aecParticles.end()) {
                                aecParticlesCells.push_back(ID);
                                sumAreaAEcParticlesCells[ID] += currentParticle->getArea();
                                aecParticles.push_back(currentParticle);
                                sumAreaAECParticles += currentParticle->getArea();
                            }

                        }
                    }
                    itP++;
                }
            }
        }

        for (int i = 0; i < 100; i++) {
            if (sumAreaAEcParticlesCells[i] > 0) {
                DEBUG_STDOUT("i have cell " << i << " with an area of " << sumAreaAEcParticlesCells[i]);
                double secretionrate = particleSecretionMoleculePerCellMin /
                                       sumAreaAEcParticlesCells[i] * time_delta;
                aecSecretionratePerGrid[i] = secretionrate;
                DEBUG_STDOUT("this cell gets an per-grid secretion rate of " +
                             std::to_string(secretionrate));
            }
        }
    }

    for (size_t i = 0; i < aecParticles.size(); i++) {
        int particleCell = aecParticlesCells[i];
        aecParticles[i]->addConcentrationChange(aecSecretionratePerGrid[particleCell]);
        if (aecParticles[i]->getId() % 100 == 0) {
        }
    }

}

std::shared_ptr<Particle> ParticleManager::createParticle(Site *site,
                                                          Coordinate3D pos,
                                                          double conc,
                                                          double dc) {
    auto particle = std::make_shared<Particle>(idHandling, site, pos, conc, dc);
    allParticles.push_back(particle);
    particleBalloonList->addCoordinateWithId(pos, idHandling);
    particleById[idHandling] = particle.get();
    noOfParticles++;
    idHandling++;

    return particle;
}

void ParticleManager::replaceParticle(Site *site, Particle *particle, std::shared_ptr<Particle> newParticle) {
    std::replace_if(allParticles.begin(),
                    allParticles.end(),
                    [particle](const auto &ptr) { return particle == ptr.get(); },
                    newParticle);
    if (newParticle == 0) noOfParticles--;
}

void ParticleManager::cleanUpAllParticles() {
    aecParticles.clear();
}

void ParticleManager::includeParticleXMLTagToc(XMLFile *xmlTags) {
    std::ostringstream ssid;
    XMLNode particlesNode = xmlTags->addChildToRootNode("Particles");
    auto itP = allParticles.begin();

    // output of particles and their neighbour-connections
    while (itP != allParticles.end()) {
        auto p = *itP;
        XMLNode particleNode = xmlTags->addChildToNode(particlesNode, "Particle");
        std::ostringstream sId, sConc, sArea, sInSite, sBoundary;
        sId << p->getId();
        Coordinate3D pos = p->getPosition();
        sConc << p->getConcentration();
        sArea << p->getArea();
        sInSite << p->getIsInSite();
        sBoundary << p->getIsAtBoundary();
        xmlTags->addDataFieldToNode(particleNode, "id", "discrete", "unsigned int", sId.str());
        xmlTags->addDataFieldToNode(particleNode, "inSite", "discrete", "unsigned int", sInSite.str());
        xmlTags->addDataFieldToNode(particleNode, "atBoundary", "discrete", "unsigned int", sBoundary.str());
        xmlTags->addDataFieldToNode(particleNode, pos);
        xmlTags->addDataFieldToNode(particleNode, "concentration", "continuous", "double", sConc.str());
        xmlTags->addDataFieldToNode(particleNode, "area", "continuous", "double", sArea.str());
        //neighbours
        XMLNode interactionPartnersNode = xmlTags->addChildToNode(particleNode, "Interactions");
        std::vector<Particle *> neighbourList = p->getParticleNeighbourList()->getNeighbours();
        std::vector<double> contactArea = p->getParticleNeighbourList()->getContactAreas();
        size_t i = 0;
        while (i < neighbourList.size()) {
            std::ostringstream sIdN, sContactN;
            sIdN << neighbourList[i]->getId();
            sContactN << contactArea[i];
            XMLNode interactionPartner = xmlTags->addChildToNode(interactionPartnersNode, "Particle");
            xmlTags->addDataFieldToNode(interactionPartner, "id", "discrete", "unsigned int", sIdN.str());
            xmlTags->addDataFieldToNode(interactionPartner, "contact", "continuous", "double", sContactN.str());
            i++;
        }
        itP++;
    }

    //output of isolines
    const int nIsolines = 5;
    double isolinesConc[nIsolines] = {24.0, 12.0, 6.0, 3, 1.5};

    auto itT = triangles.begin();
    while (itT != triangles.end()) {
        TRIANGLE3D t3d = *itT;
        auto p1 = allParticles[t3d.neighbourIds[0]];
        auto p2 = allParticles[t3d.neighbourIds[1]];
        auto p3 = allParticles[t3d.neighbourIds[2]];

        double c1, c2, c3;
        c1 = p1->getConcentration();
        c2 = p2->getConcentration();
        c3 = p3->getConcentration();

        double minConc = std::min(std::min(c1, c2), c3);
        double maxConc = std::max(std::max(c1, c2), c3);

        Coordinate3D pointIso1 = Coordinate3D(), pointIso2 = Coordinate3D();
        for (int i = 0; i < nIsolines; i++) {

            if (minConc <= isolinesConc[i] && maxConc >= isolinesConc[i]) {

                //find the positions of isoline vertices on the triangle
                double relDist12 = (isolinesConc[i] - c1) / (c2 - c1);
                double relDist23 = (isolinesConc[i] - c2) / (c3 - c2);
                double relDist31 = (isolinesConc[i] - c3) / (c1 - c3);

                if (relDist12 < 1 && relDist12 > 0) {
                    Coordinate3D connect(p2->getPosition() - p1->getPosition());
                    connect *= relDist12;
                    pointIso1 = p1->getPosition();
                    pointIso1 += connect;
                }

                if (relDist23 < 1 && relDist23 > 0) {
                    Coordinate3D connect(p3->getPosition() - p2->getPosition());
                    connect *= relDist23;

                    if (pointIso1.x == 0 && pointIso1.y == 0 && pointIso1.z == 0) {
                        pointIso1 = p2->getPosition();
                        pointIso1 += connect;
                    } else {
                        pointIso2 = p2->getPosition();
                        pointIso2 += connect;
                    }
                }

                if (relDist31 < 1 && relDist31 > 0) {
                    Coordinate3D connect(p1->getPosition() - p3->getPosition());
                    connect *= relDist31;
                    pointIso2 = p3->getPosition();
                    pointIso2 += connect;
                }
            }

            //add to xml file
            bool isZeroIso1OrIso2 = (pointIso1.x == 0 && pointIso1.y == 0 && pointIso1.z == 0) || (pointIso2.x == 0 && pointIso2.y == 0 && pointIso2.z == 0);
            if (!isZeroIso1OrIso2) {
                XMLNode isolineNode = xmlTags->addChildToNode(particlesNode, "Isoline");
                xmlTags->addDataFieldToNode(isolineNode, pointIso1, "start");
                xmlTags->addDataFieldToNode(isolineNode, pointIso2, "end");
                std::ostringstream sInt;
                sInt << i;
                xmlTags->addDataFieldToNode(isolineNode, "strength", "continuous", "double", sInt.str());
            }
            pointIso1 = Coordinate3D();
            pointIso2 = Coordinate3D();
        }
        itT++;
    }
}

void ParticleManager::diffusionPSE(double timestep, double time_delta) {
    for (auto itP = allParticles.begin(); itP != allParticles.end(); itP++) {
        (*itP)->diffusePSE(timestep);
    }
    double maxChange = 0, change;
    for (auto itP = allParticles.begin(); itP != allParticles.end(); itP++) {
        change = (*itP)->applyConcentrationChange(time_delta);
        if (change > maxChange) {
            maxChange = change;
        }
    }
}

bool ParticleManager::steadyStateReached(double current_time) {
    bool stStReached = false;
    if (dc > 500) { //below 500, reaching a steady state takes too much time
        // all values below were found experimentally
        if (dc == 600) {
            stStReached =
                    (current_time - site->getAgentManager()->getLastConidiaChange() > 35.0);
            allowHigherDT = stStReached;
        } else if (dc == 2000) {
            stStReached =
                    (current_time - site->getAgentManager()->getLastConidiaChange() > 13.0);
            allowHigherDT = stStReached;
        } else if (dc == 6000) {
            stStReached =
                    (current_time - site->getAgentManager()->getLastConidiaChange() > 5.0);
            allowHigherDT = stStReached;
        } else if (dc == 15000) {
            stStReached =
                    (current_time - site->getAgentManager()->getLastConidiaChange() > 2.2);
            allowHigherDT = stStReached;
        } else if (dc == 35000) {
            stStReached =
                    (current_time - site->getAgentManager()->getLastConidiaChange() > 1.0);
            allowHigherDT = stStReached;
        } else {
            double limit = 30000.0/dc; // 30000.0/DC proved to be a valid value until steady state is reached for high DC
            stStReached = (current_time - site->getAgentManager()->getLastConidiaChange() > limit);
            allowHigherDT = stStReached;
        }
    }
    return stStReached;
}

double ParticleManager::getGradient(const Coordinate3D &position) {
    double gradient{};
    std::vector<unsigned int> closest3ParticleIdxs;
    this->getParticleBalloonList()->getClosestObjectIndices(position, closest3ParticleIdxs, 3);
    auto p1 = allParticles[closest3ParticleIdxs[0]];
    auto p2 = allParticles[closest3ParticleIdxs[1]];
    auto p3 = allParticles[closest3ParticleIdxs[2]];
    double linearlyInterpolatedValue = Algorithms::interpolateBilinearOnTriangle(position, p1->getPosition(),
                                                                                 p1->getGradientStrength(),
                                                                                 p2->getPosition(),
                                                                                 p2->getGradientStrength(),
                                                                                 p3->getPosition(),
                                                                                 p3->getGradientStrength());
    gradient = linearlyInterpolatedValue;
    return gradient;
}

void ParticleManager::extractTriangles() {
    auto itP = allParticles.begin();
    while (itP != allParticles.end()) {
        auto p = *itP;
        int idP = p->getId();
        size_t i = 0, j = 0;
        std::vector<Particle *> neighbourList = p->getParticleNeighbourList()->getNeighbours();
        while (i < neighbourList.size()) {
            j = i + 1;
            int idI = neighbourList[i]->getId();;
            //check if neighbours (i,j) are neighbours of each other
            while (j < neighbourList.size()) {
                int idJ = neighbourList[j]->getId();
                if (allParticles[idI]->getParticleNeighbourList()->existsInList(neighbourList[j])) {
                    addTriangle(idP, idI, idJ);
                }
                j++;
            }
            i++;
        }
        itP++;
    }
}

void ParticleManager::addTriangle(unsigned int id1, unsigned int id2, unsigned int id3) {
    std::vector<unsigned int> sorted;
    sorted.push_back(id1);
    sorted.push_back(id2);
    sorted.push_back(id3);
    std::sort(sorted.begin(), sorted.end());
    bool toBeAdded = true;

    //check if this triangle is already in the list of triangles
    auto itT = triangles.begin();
    while (itT != triangles.end()) {
        TRIANGLE3D t3d = *itT;
        if (t3d.neighbourIds[0] == sorted[0]) {
            if (t3d.neighbourIds[1] == sorted[1]) {
                if (t3d.neighbourIds[2] == sorted[2]) {
                    toBeAdded = false;
                }
            }
        }

        itT++;
    }

    if (toBeAdded) {
        TRIANGLE3D newTriangle{};
        newTriangle.outside = TRUE;
        newTriangle.neighbourIds[0] = sorted[0];
        newTriangle.neighbourIds[1] = sorted[1];
        newTriangle.neighbourIds[2] = sorted[2];
        if (allParticles[newTriangle.neighbourIds[0]]->getIsInSite()) {
            newTriangle.outside = FALSE;
        } else {
            if (allParticles[newTriangle.neighbourIds[1]]->getIsInSite()) {
                newTriangle.outside = FALSE;
            } else {
                if (allParticles[newTriangle.neighbourIds[2]]->getIsInSite()) {
                    newTriangle.outside = FALSE;
                }
            }
        }
        triangles.push_back(newTriangle);
    }
}

double ParticleManager::getSumChemokine() {
    auto itP = allParticles.begin();
    double sumOfChemokine = 0.0;
    int counter = 0;
    while (itP != allParticles.end()) {
        auto p = *itP;
        double conc_p = p->getConcentration();
        sumOfChemokine += conc_p;
        itP++;
    }
    return sumOfChemokine;
}

int ParticleManager::get_closest_AEC_ID(Coordinate3D position, int type) {    //position of particle and cell type

    if (type == 1) {
        int id = -1;
        double mindist = 1000000;
        for (size_t i = 0; i < alvEpithTypeOne.size(); i++) {
            double dist = alvEpithTypeOne.at(i).calculateEuclidianDistance(abm::util::toSphericCoordinates(position));
            if (dist < mindist) {
                mindist = dist;
                id = i;
            }
        }
        return id;
    }
    if (type == 2) {
        int id = -1;
        double mindist = 1000000;
        for (size_t i = 0; i < alvEpithTypeTwo.size(); i++) {
            double dist = alvEpithTypeTwo.at(i).calculateEuclidianDistance(abm::util::toSphericCoordinates(position));
            if (dist < mindist) {
                mindist = dist;
                id = i;
            }
        }
        return id;
    }
    return 0;
}

void ParticleManager::setAECCells(std::vector<SphericCoordinate3D> AECT1, std::vector<SphericCoordinate3D> AECT2) {
    alvEpithTypeOne = AECT1;
    alvEpithTypeTwo = AECT2;
    DEBUG_STDOUT("particle Manager got information about " << std::to_string(alvEpithTypeOne.size()) << " AECT1 cells");
    DEBUG_STDOUT("particle Manager got information about " << std::to_string(alvEpithTypeTwo.size()) << " AECT2 cells");
}
