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

#include "SphericalMorphology.h"
#include "SphereRepresentation.h"
#include "simulation/Site.h"

SphericalMorphology::SphericalMorphology() : MorphologyElement() {
}

SphericalMorphology::SphericalMorphology(Morphology *refMorphology, std::shared_ptr<Coordinate3D> position,
                                         double radius, std::string description) : MorphologyElement(refMorphology,
                                                                                                     description) {
    center_of_mass_ = position;
    this->radius = radius;
    this->description_ = description;

    generateSphereRepresentation();
}

void SphericalMorphology::generateSphereRepresentation() {
    auto sR = std::make_unique<SphereRepresentation>(center_of_mass_, radius, this, description_);
    sphere_rep_.push_back(std::move(sR));
}

std::string SphericalMorphology::generatePovObject() {
    std::ostringstream ss;
    ss << "sphere {" << '\n';
    ss << "\t" << center_of_mass_->x << "," << center_of_mass_->y << "," << center_of_mass_->z << "," << '\n';
    ss << "\t" << getRadius() << '\n';
    ss << "\t" << "pigment { color " << morphology_this_belongs_to_->getColorRGB()->printPovColorRGBT() << " }" << '\n';
    ss << "\t}" << '\n';

    return ss.str();
}

void SphericalMorphology::includeMorphologyXMLOutput(XMLFile *xmlFile, XMLNode *node) {
    XMLNode sphericalMorph = xmlFile->addChildToNode(*node, "SphericalMorphology");
    std::ostringstream radiusStr;
    radiusStr << getRadius();
    xmlFile->addDataFieldToNode(sphericalMorph, "radius", "continuous", "double", radiusStr.str());
}

void SphericalMorphology::includeAssociatedMorphologyXMLOutput(XMLFile *xmlFile, XMLNode *morphNode) {
    XMLNode sphericalMorph = xmlFile->addChildToNode(*morphNode, "SphericalMorphology");
    std::ostringstream radiusStr;
    radiusStr << getRadius();
    xmlFile->addDataFieldToNode(sphericalMorph, "radius", "continuous", "double", radiusStr.str());
    xmlFile->addDataFieldToNode(sphericalMorph, *center_of_mass_);
}