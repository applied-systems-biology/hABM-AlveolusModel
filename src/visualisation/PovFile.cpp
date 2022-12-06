//  Copyright by Christoph Saffer, Paul Rudolph, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#include "visualisation/PovFile.h"
#include "utils/io_util.h"
#include "simulation/Site.h"
#include "visualisation/PovRayObject.h"
#include "utils/macros.h"

PovFile::PovFile(const std::string &filebody) {
    filebody_pov_ = filebody;
}

void PovFile::setGlobalPart(std::string global) {
    global_part_ = global;
}

void PovFile::setCameraPart(std::string cam) {
    camera_part_ = cam;
}

void PovFile::setBackgroundPart(std::string bg) {
    background_part_ = bg;
}

void PovFile::setLightPart(std::string light) {
    light_part_ = light;
}

void PovFile::addLightPart(const std::string &light) {
    light_part_ += light;
}

void PovFile::setBorderPart(std::string border) {
    border_part_ = border;
}

void PovFile::addPovObject(const std::string &povObj) {
    pov_objects_.emplace_back(povObj);
}

void PovFile::setDimensions(const std::string px_width, const std::string px_height) {
    px_width_ = px_width, px_height_ = px_height;
}

void PovFile::doPovProcess(double current_time, std::string zpos) {
    std::ofstream output(filebody_pov_ + zpos + ".pov");
    output << global_part_ << '\n';
    output << camera_part_ << '\n';
    output << light_part_ << '\n';
    output << border_part_ << '\n';
    output << background_part_ << '\n';
    for (const auto &obj: pov_objects_) {
        output << obj << '\n';
    }
    bool include_time = false;
    output.close();
    int left_padding = static_cast<int>(stoi(px_width_) / 6);
    int below_padding = 50;
    int font_size = static_cast<int>(stoi(px_height_) / 30);
    std::string font_color = "black";
    int X = std::stoi(px_width_) - left_padding;
    int Y = std::stoi(px_height_) - below_padding;
    char timeOfFrame[20];
    sprintf(timeOfFrame, "%.3f", current_time);
    std::ostringstream command;
    command << "povray -GA -d +H" << px_height_ << " +W" << px_width_ << " +I" << filebody_pov_ << zpos << ".pov +O"
            << filebody_pov_ << zpos << ".png >/dev/null 2>&1 ";
    if (include_time) {
        command << "&& convert -pointsize " << font_size << " -fill " << font_color << " -draw \"text " << X << " " << Y
                << " '" << timeOfFrame
                << " min'\" " << filebody_pov_ << zpos << ".png " << filebody_pov_ << ".png ";
    }
    command << "&& rm -f " << filebody_pov_ << zpos << ".pov";

    abm::util::executeShellCommand(command.str(), true);
}

void PovFile::transcribeSite(const Site &site) {
    auto file = std::make_unique<XMLFile>();
    site.includeSiteXMLTagToc(file.get());
    auto root = file->getRootNode();
    int nodes = root.nChildNode();
    for (int i = 0; i < nodes; i++) {

        XMLNode child = root.getChildNode(i);
        std::string nameOfNode_site = child.getName();
        if (nameOfNode_site == "Site") {
            auto specificSiteNode = child.getChildNode(0);
            std::string nameOfNode = specificSiteNode.getName();

            if (nameOfNode == "SphereSite") {
                Coordinate3D position = InputConfiguration::getCoordinateDataFieldValueByName(&specificSiteNode,
                                                                                              "position");
                double radius = InputConfiguration::getDoubleDataFieldValueByName(&specificSiteNode, "radius");
                ColorRGB crgb = ColorRGB(1, 0, 0, 0.7);
                ColorRGB crgbTori = ColorRGB(0.5, 0.5, 0.5, 0.5);
                addPovObject(PovRayObject::getSphere(position, radius, crgb));
                addPovObject(PovRayObject::getTorus(position, 1.5, radius, Coordinate3D{0, -90, 0}, &crgbTori));
                addPovObject(PovRayObject::getTorus(position, 1.5, radius, Coordinate3D{-90, 0, 0}, &crgbTori));
                addPovObject(PovRayObject::getTorus(position, 1.5, radius, Coordinate3D{0, 0, -90}, &crgbTori));
            }

            if (nameOfNode == "AlveoleSite") {
                Coordinate3D position = InputConfiguration::getCoordinateDataFieldValueByName(&specificSiteNode,
                                                                                              "position");
                double radius = InputConfiguration::getDoubleDataFieldValueByName(&specificSiteNode, "radius");
                double thetaLowerBound = InputConfiguration::getDoubleDataFieldValueByName(&specificSiteNode,
                                                                                           "theta-lower-bound");
                double thickness = InputConfiguration::getDoubleDataFieldValueByName(&specificSiteNode, "thickness");
                ColorRGB crgbE1 = ColorRGB(1, 1, 0, 0.8);
                ColorRGB crgbE2 = ColorRGB(0, 1, 1, 0.8);
                ColorRGB crgbPoK = ColorRGB(0, 0, 0, 0.5);
                ColorRGB crgbTori = ColorRGB(0.5, 0.5, 0.5, 0.5);

                //big box induced by the lower bound of theta:
                Coordinate3D rightUpBox = Coordinate3D{radius + thickness, radius + thickness,
                                                       radius * cos(thetaLowerBound)};
                Coordinate3D lowDownBox = Coordinate3D{-(radius + thickness), -(radius + thickness),
                                                       -(radius + thickness)};
                std::string povBoxAlveolusBoundaries = PovRayObject::getBox(lowDownBox, rightUpBox, crgbE1);

                int numberOfNodes = specificSiteNode.nChildNode();
                for (int i = 0; i < numberOfNodes; i++) {
                    XMLNode curNode = specificSiteNode.getChildNode(i);
                    std::string nameOfNode_s = curNode.getName();

                    if (nameOfNode_s == "AET1") {

                        double radiusCell = InputConfiguration::getDoubleDataFieldValueByName(&curNode, "radius");
                        Coordinate3D posCell = InputConfiguration::getCoordinateDataFieldValueByName(&curNode,
                                                                                                     "position");

                        std::string cellSphere = PovRayObject::getSphere(posCell, radiusCell, crgbE1);
                        std::string outerBoundary = PovRayObject::getSphere(position, radius + thickness, crgbE1);
                        std::string innerBoundary = PovRayObject::getSphere(position, radius, crgbE1);
                        std::string boundarySite = PovRayObject::getDifference(outerBoundary, innerBoundary);
                        std::string cellPlate = PovRayObject::getIntersection(boundarySite, cellSphere);
                        std::vector<std::string> boxes; //all neighbouring cells including the big box induced by the lower bound of theta
                        std::vector<std::string> boxEC2;

                        boxes.push_back(povBoxAlveolusBoundaries);

                        int numberOfSubNodes = curNode.nChildNode();
                        for (int j = 0; j < numberOfSubNodes; j++) {

                            XMLNode curSubNode = curNode.getChildNode(j);
                            std::string nameOfSubNode = curSubNode.getName();

                            if (nameOfSubNode == "AETNC1") {
                                SphericCoordinate3D nConnect = InputConfiguration::getSphericVectorDataFieldByName(
                                        &curSubNode, "position");
                                Coordinate3D rightUp = Coordinate3D{radiusCell, radiusCell, nConnect.r / 2.0};
                                Coordinate3D lowDown = Coordinate3D{-radiusCell, -radiusCell, -radiusCell};
                                boxes.push_back(PovRayObject::getBox(lowDown, rightUp, crgbE1, nConnect, posCell));
                            }

                            if (nameOfSubNode == "AETNC2") {
                                Coordinate3D posBox = InputConfiguration::getCoordinateDataFieldValueByName(&curSubNode,
                                                                                                            "position");
                                double lengthCell = InputConfiguration::getDoubleDataFieldValueByName(&curSubNode,
                                                                                                      "length");
                                SphericCoordinate3D posCellSpheric = abm::util::toSphericCoordinates(posBox);
                                SphericCoordinate3D
                                        posCellSphericVec = SphericCoordinate3D{posCellSpheric.r, posCellSpheric.theta,
                                                                                posCellSpheric.phi};
                                Coordinate3D rightUp = Coordinate3D{lengthCell / 2.0, lengthCell / 2.0,
                                                                    lengthCell / 2.0};
                                Coordinate3D lowDown = Coordinate3D{-lengthCell / 2.0, -lengthCell / 2.0,
                                                                    -lengthCell / 2.0};
                                boxEC2.push_back(
                                        PovRayObject::getBox(lowDown, rightUp, crgbE1, posCellSphericVec, posBox));
                            }

                            if (nameOfSubNode == "PoKNC") {
                                Coordinate3D posPoK = InputConfiguration::getCoordinateDataFieldValueByName(&curSubNode,
                                                                                                            "position");
                                double radiusPoK = InputConfiguration::getDoubleDataFieldValueByName(&curSubNode,
                                                                                                     "radius");
                                SphericCoordinate3D sphericPosPok = abm::util::toSphericCoordinates(posPoK);
                                sphericPosPok.r += thickness;
                                Coordinate3D posPoK2 = abm::util::toCartesianCoordinates(sphericPosPok);
                                sphericPosPok.r -= sphericPosPok.r - (thickness + 1);
                                posPoK = abm::util::toCartesianCoordinates(sphericPosPok);
                                //cout << "[TTPT] position PoK " << posPoK.printCoordinate() << " radius " << radiusPoK << endl;

                                std::string poreSphere = PovRayObject::getCylinder(posPoK, posPoK2, radiusPoK, crgbPoK);

                                boxEC2.push_back(poreSphere);

                            }
                        }
                        std::string voronoiCell = PovRayObject::getIntersection(cellPlate, boxes);
                        std::string voronoiCellDiffBoxes = PovRayObject::getDifference(voronoiCell, boxEC2);
                        //cout << "[TTPT] " << voronoiCell << endl;
                        addPovObject(voronoiCellDiffBoxes);
                    }

                    if (nameOfNode_s == "AET2") {
                        double lengthCell = InputConfiguration::getDoubleDataFieldValueByName(&curNode, "length");
                        Coordinate3D posBox = InputConfiguration::getCoordinateDataFieldValueByName(&curNode,
                                                                                                    "position");
                        SphericCoordinate3D posCellSpheric = abm::util::toSphericCoordinates(posBox);
                        SphericCoordinate3D posCellSphericVec = SphericCoordinate3D{posCellSpheric.r,
                                                                                    posCellSpheric.theta,
                                                                                    posCellSpheric.phi};
                        Coordinate3D rightUp = Coordinate3D{lengthCell / 2.0, lengthCell / 2.0, lengthCell / 2.0};
                        Coordinate3D lowDown = Coordinate3D{-lengthCell / 2.0, -lengthCell / 2.0, -lengthCell / 2.0};
                        std::string cellBox = PovRayObject::getBox(lowDown, rightUp, crgbE2, posCellSphericVec, posBox);
                        std::string outerBoundary = PovRayObject::getSphere(position, radius + thickness, crgbE2);
                        std::string innerBoundary = PovRayObject::getSphere(position, radius, crgbE2);
                        std::string boundarySite = PovRayObject::getDifference(outerBoundary, innerBoundary);
                        std::string boundarySite2 = PovRayObject::getIntersection(boundarySite,
                                                                                  povBoxAlveolusBoundaries);
                        std::string boxPlate = PovRayObject::getIntersection(boundarySite2, cellBox);
                        addPovObject(boxPlate);
                    }
                    if (nameOfNode_s == "PoK") {
                        double radiusPoK = InputConfiguration::getDoubleDataFieldValueByName(&curNode, "radius");
                        Coordinate3D posPoK = InputConfiguration::getCoordinateDataFieldValueByName(&curNode,
                                                                                                    "position");
                        SphericCoordinate3D sphericPosPok = abm::util::toSphericCoordinates(posPoK);
                        sphericPosPok.r += thickness;
                        Coordinate3D posPoK2 = abm::util::toCartesianCoordinates(sphericPosPok);
                        sphericPosPok.r -= (thickness + 1);
                        posPoK = abm::util::toCartesianCoordinates(sphericPosPok);
                        std::string poreSphere = PovRayObject::getCylinder(posPoK, posPoK2, radiusPoK, crgbPoK);
                        std::string outerBoundary = PovRayObject::getSphere(position, radius + thickness, crgbPoK);
                        std::string innerBoundary = PovRayObject::getSphere(position, radius, crgbPoK);
                        std::string boundarySite = PovRayObject::getDifference(outerBoundary, innerBoundary);
                        std::string platePoreOfKohn = PovRayObject::getIntersection(boundarySite, poreSphere);
                        addPovObject(platePoreOfKohn);
                    }
                }
            }
        }
    }
}

void PovFile::transcribeParticles(const Site &site) {
    constexpr auto globalFade = 1.0; //1.0 - particles visible, 0.0 - site visible
    auto file = std::make_unique<XMLFile>();
    site.getParticleManager()->includeParticleXMLTagToc(file.get());
    auto root = file->getRootNode();
    int nodes = root.nChildNode();
    for (int i = 0; i < nodes; i++) {

        XMLNode curNode = root.getChildNode(i);
        std::string nameOfNode = curNode.getName();

        if (nameOfNode == "Particles") {
            int n = curNode.nChildNode();
            double xVal = 0, yVal = 0, zVal = 0;
            std::vector<Coordinate3D> positions, voronois;

            Coordinate3D startIso, endIso;
            double strength;

            for (int j = 0; j < n; j++) {
                XMLNode particleNode = curNode.getChildNode(j);
                nameOfNode = particleNode.getName();

                if (nameOfNode == "Particle") {
                    Coordinate3D position = InputConfiguration::getCoordinateDataFieldValueByName(&particleNode,
                                                                                                  "position");
                    positions.push_back(position);
                }

                if (nameOfNode == "Voronoi") {
                    Coordinate3D position = InputConfiguration::getCoordinateDataFieldValueByName(&particleNode,
                                                                                                  "position");
                    voronois.push_back(position);
                }

                if (nameOfNode == "Isoline") {

                    startIso = InputConfiguration::getCoordinateDataFieldValueByName(&particleNode, "start");
                    endIso = InputConfiguration::getCoordinateDataFieldValueByName(&particleNode, "end");
                    strength = InputConfiguration::getDoubleDataFieldValueByName(&particleNode, "strength");

                    ColorRGB crgb = ColorRGB(1.0, 1.0, 1.0, 0.0);
                    double radius = 1.5 / (strength + 1);
                    addPovObject(PovRayObject::getCylinder(startIso, endIso, radius, crgb));
                    addPovObject(PovRayObject::getSphere(startIso, radius, crgb));
                    addPovObject(PovRayObject::getSphere(endIso, radius, crgb));

                }
            }

            for (int j = 0; j < n; j++) {
                XMLNode particleNode = curNode.getChildNode(j);
                nameOfNode = particleNode.getName();

                if (nameOfNode == "Particle") {
                    Coordinate3D position = positions[j];
                    ColorRGB crgb;
                    double conc = InputConfiguration::getUIntDataFieldValueByName(&particleNode, "concentration");
                    if (InputConfiguration::getUIntDataFieldValueByName(&particleNode, "inSite") == 0) {

                        if (InputConfiguration::getUIntDataFieldValueByName(&particleNode, "atBoundary") == 0) {
                            continue; //skip this particle
                        } else {
                            crgb = ColorRGB(0.0, 0.0, 1.0, 0.5 + ((1.0 - globalFade) * 0.5)); //blue
                        }
                    } else {
                        crgb = ColorRGB(conc / 20.0f, 0.0, conc / 20.0f, 0.5 + ((1.0 - globalFade) * 0.5)); //orange
                    }
                    double siteRadius = site.getRadius();
                    double radius = 0.75;
                    SphericCoordinate3D sc3d = abm::util::toSphericCoordinates(position);
                    sc3d.r = sc3d.r;
                    position = abm::util::toCartesianCoordinates(sc3d);
                    sc3d.r -= 0.1;
                    Coordinate3D iPostion = abm::util::toCartesianCoordinates(sc3d);
                    if (InputConfiguration::getUIntDataFieldValueByName(&particleNode, "inSite") != 0) {
                        XMLNode interactionsNode = particleNode.getChildNode("Interactions");
                        for (int k = 0; k < interactionsNode.nChildNode(); k++) {
                            XMLNode childInteractionNode = interactionsNode.getChildNode(k);
                            std::string nameOfChild = childInteractionNode.getName();
                            if (nameOfChild == "Particle") {
                                unsigned int idN = InputConfiguration::getUIntDataFieldValueByName(
                                        &childInteractionNode, "id");
                                Coordinate3D positionN = positions[idN];
                                double radiusCyl = 0.13;
                            }
                        }
                    }
                }
                unsigned int particles = positions.size();
                if (nameOfNode == "Voronoi") {
                    Coordinate3D position = voronois[j - particles];
                    ColorRGB crgb = ColorRGB(0.0, 1.0, 0.0, 0.5 + ((1.0 - globalFade) * 0.5));
                    double radius = 0.2;
                    XMLNode interactionsNode = particleNode.getChildNode("Interactions");
                    for (int k = 0; k < interactionsNode.nChildNode(); k++) {
                        XMLNode childInteractionNode = interactionsNode.getChildNode(k);
                        std::string nameOfChild = childInteractionNode.getName();
                        if (nameOfChild == "Voronoi") {
                            unsigned int idN = InputConfiguration::getUIntDataFieldValueByName(&childInteractionNode,
                                                                                               "id");
                            Coordinate3D positionN = voronois[idN];
                            ColorRGB crgbCyl = ColorRGB(0.0, 1, 0.0, 0.5 + ((1.0 - globalFade) * 0.5));
                            double radiusCyl = 0.2;
                        }
                    }
                }
            }
        }
    }
}

void PovFile::transcribeAgents(const Site &site) {
    // set background color
    addPovObject("  background { color rgb <0.9, 0.9, 0.9> }");
    auto file = std::make_unique<XMLFile>();
    site.getAgentManager()->getAgentXMLTocTags(file.get(), false);
    auto root = file->getRootNode();
    double site_radius = site.getRadius();
    int nodes = root.nChildNode();
    for (int i = 0; i < nodes; i++) {

        XMLNode curAgentNode = root.getChildNode(i);
        std::string nameOfNode = curAgentNode.getName();

        if (nameOfNode == "Agent") {
            XMLNode curCellNode = curAgentNode.getChildNode("Cell");
            XMLNode curMorphology = curCellNode.getChildNode("Morphology");
            ColorRGB crgb = InputConfiguration::getColorDataFieldByName(&curMorphology, "color");
            std::string agentname =
                    InputConfiguration::getStringDataFieldByName(&curAgentNode, "typename");
            //            cout << "transcribung agent: " << agentname << endl;

            // loop over all morphology types involved in this cell
            int morphNodes = curMorphology.nChildNode();
            for (int k = 0; k < morphNodes; k++) {
                XMLNode curMorphPartNode = curMorphology.getChildNode(k);
                Coordinate3D position = InputConfiguration::getCoordinateDataFieldValueByName(&curMorphPartNode,
                                                                                              "position");
                double radius = InputConfiguration::getDoubleDataFieldValueByName(&curMorphPartNode, "radius");
                std::string nameOfMorphNode = curMorphPartNode.getName();

                if (nameOfMorphNode == "SphericalMorphology") {
                    if (site_radius > 0) { // Since a radius greater 0 exists, we are in AlveoleSite
                        if (agentname == "AspergillusFumigatus") {
                            if (abm::util::toSphericCoordinates(position).r > site_radius * 1.025) {
                                crgb = {1, 1, 0.0};
                            } else {
                                crgb = InputConfiguration::getColorDataFieldByName(&curMorphology, "color");
                                addPovObject(PovRayObject::getSphere(position, radius, crgb));
                            }
                        } else {
                            double radius = InputConfiguration::getDoubleDataFieldValueByName(&curMorphPartNode,"radius");
                            std::string AMsphere = PovRayObject::getSphere(position, radius, crgb);
                            Coordinate3D center = Coordinate3D{0, 0, 0};
                            std::string alveole = PovRayObject::getSphere(center, site_radius, crgb);
                            std::string smallalveole = PovRayObject::getSphere(center, site_radius - 2, crgb);
                            std::string innenAM = PovRayObject::getIntersection(AMsphere, alveole);
                            addPovObject(PovRayObject::getDifference(innenAM, smallalveole));
                        }
                    }
                }

                if (nameOfMorphNode == "CylindricMorphology") {
                    double radius = InputConfiguration::getDoubleDataFieldValueByName(&curMorphPartNode, "radius");
                    Coordinate3D start = InputConfiguration::getCoordinateDataFieldValueByName(&curMorphPartNode,
                                                                                               "start");
                    Coordinate3D end = InputConfiguration::getCoordinateDataFieldValueByName(&curMorphPartNode, "end");

                    if (InputConfiguration::getUIntDataFieldValueByName(&curMorphPartNode, "smoothingOn") == 1) {
                        addPovObject(PovRayObject::getCylinder(start, end, radius, crgb));
                        addPovObject(PovRayObject::getSphere(start, radius, crgb));
                        addPovObject(PovRayObject::getSphere(end, radius, crgb));
                    } else {
                        addPovObject(PovRayObject::getCylinder(start, end, radius, crgb));
                    }
                }
            }
        }
    }
}

void PovFile::transcribeAgents(const Site &site, double z_start, double z_end) {

    auto file = std::make_unique<XMLFile>();
    site.getAgentManager()->getAgentXMLTocTags(file.get(), false);
    auto root = file->getRootNode();
    ColorRGB boxColor = ColorRGB(0, 1, 0, 0.85); //green box
    Coordinate3D leftDown = Coordinate3D{-1000, -1000, z_start};
    Coordinate3D rightUp = Coordinate3D{1000, 1000, z_end};
    std::string cutBox = PovRayObject::getBox(leftDown, rightUp, boxColor);

    int nodes = root.nChildNode();
    for (int i = 0; i < nodes; i++) {

        XMLNode curAgentNode = root.getChildNode(i);
        std::string nameOfNode = curAgentNode.getName();

        if (nameOfNode == "Agent") {

            Coordinate3D position = InputConfiguration::getCoordinateDataFieldValueByName(&curAgentNode, "position");

            XMLNode curCellNode = curAgentNode.getChildNode("Cell");
            XMLNode curMorphology = curCellNode.getChildNode("Morphology");
            ColorRGB crgb = InputConfiguration::getColorDataFieldByName(&curMorphology, "color");

            //loop over all morphology types involved in this cell
            int morphNodes = curMorphology.nChildNode();
            for (int k = 0; k < morphNodes; k++) {
                XMLNode curMorphPartNode = curMorphology.getChildNode(k);
                Coordinate3D position = InputConfiguration::getCoordinateDataFieldValueByName(&curMorphPartNode,
                                                                                              "position");
                std::string name_of_morph_node = curMorphPartNode.getName();

                if (name_of_morph_node == "SphericalMorphology") {
                    double radius = InputConfiguration::getDoubleDataFieldValueByName(&curMorphPartNode, "radius");
                    double zPos = position.z;
                    if (((zPos > z_start) && (zPos < z_end))) {
                        std::string sphere = PovRayObject::getSphere(position, radius, crgb);
                        addPovObject(PovRayObject::getIntersection(sphere, cutBox));
                    } else {
                        if (((zPos + radius > z_start) && (zPos + radius < z_end)) ||
                            ((zPos - radius > z_start) && (zPos - radius < z_end))) {
                            std::string sphere = PovRayObject::getSphere(position, radius, crgb);
                            addPovObject(PovRayObject::getIntersection(sphere, cutBox));
                        }
                    }
                }
                if (name_of_morph_node == "CylindricMorphology") {
                    double radius = InputConfiguration::getDoubleDataFieldValueByName(&curMorphPartNode, "radius");
                    Coordinate3D start = InputConfiguration::getCoordinateDataFieldValueByName(&curMorphPartNode,
                                                                                               "start");
                    Coordinate3D end = InputConfiguration::getCoordinateDataFieldValueByName(&curMorphPartNode, "end");

                    if (InputConfiguration::getUIntDataFieldValueByName(&curMorphPartNode, "smoothingOn") == 1) {
                        addPovObject(PovRayObject::getCylinder(start, end, radius, crgb));
                        addPovObject(PovRayObject::getSphere(start, radius, crgb));
                        addPovObject(PovRayObject::getSphere(end, radius, crgb));
                    } else {
                        addPovObject(PovRayObject::getCylinder(start, end, radius, crgb));
                    }
                }
            }
        }
    }
}