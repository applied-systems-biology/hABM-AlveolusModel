//  Copyright by Christoph Saffer, Paul Rudolph, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#include <list>

#include "simulation/Morphology.h"
#include "simulation/Cell.h"
#include "simulation/morphology/SphereRepresentation.h"

#ifndef M_PI
#define M_PI    3.14159265358979323846f
#endif


Morphology::Morphology(const std::string &color, Cell *cell) {
    color_rgb_ = std::make_unique<ColorRGB>(cell->getSite()->getRandomGenerator(), color);
    this->cell_this_belongs_to_ = cell;
}

std::string Morphology::generatePovObject() {
    return "";
}

double Morphology::getMaximumDistanceFromCenter() {
    return 99999;
}

void Morphology::includeMorphologyXMLOutput(XMLFile *xmlFile, XMLNode *node) {
    XMLNode morphologyNode = xmlFile->addChildToNode(*node, "Morphology");
    xmlFile->addDataFieldToNode(morphologyNode, *color_rgb_);

    auto it = morphologyElements.begin();
    while (it != morphologyElements.end()) {
        (*it)->includeAssociatedMorphologyXMLOutput(xmlFile, &morphologyNode);
        it++;
    }
}

std::string Morphology::getTypeName() {
    return "Morphology";
}

ColorRGB *Morphology::getColorRGB() {
    return color_rgb_.get();
}

void Morphology::setColorRGB(std::unique_ptr<ColorRGB> col) {
    color_rgb_ = std::move(col);
}

void Morphology::appendAssociatedCellpart(std::unique_ptr<MorphologyElement> morphElement) {
    morphologyElements.emplace_back(std::move(morphElement));
}

std::vector<SphereRepresentation *> Morphology::getAllSpheresOfThis() {
    std::vector<SphereRepresentation *> allSpheresOfThisCell;

    auto it = morphologyElements.begin();

    while (it != morphologyElements.end()) {
        std::vector<std::shared_ptr<SphereRepresentation>> sphereReps = (*it)->getSphereRepresentation();
        for (const auto &ptr: sphereReps) {
            allSpheresOfThisCell.emplace_back(ptr.get());
        }
        it++;
    }

    return allSpheresOfThisCell;
}

SphereRepresentation *Morphology::getBasicSphereOfThis() {
    std::vector<SphereRepresentation *> allSpheresOfThisCell = getAllSpheresOfThis();
    auto it = allSpheresOfThisCell.begin();
    SphereRepresentation *basicSphere = 0;
    while (it != allSpheresOfThisCell.end()) {
        if ((*it)->getMorphologyElementThisBelongsTo()->getDescription() == "basic") {
            basicSphere = *it;
            break;
        }
        it++;
    }

    return basicSphere;
}

double Morphology::getVolume() {
    double volume = 0;
    std::vector<SphereRepresentation *> spheres = getAllSpheresOfThis();
    for (auto &sphere : spheres) {
        volume += (4 / 3) * M_PI * pow(sphere->getRadius(), 3);
    }
    return volume;
}