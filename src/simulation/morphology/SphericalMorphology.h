//  Copyright by Christoph Saffer, Paul Rudolph, Sandra Timme, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#ifndef SPHERESURFACE_H
#define    SPHERESURFACE_H

#include "MorphologyElement.h"
#include "SphereRepresentation.h"

class SphericalMorphology : public MorphologyElement {
public:
  // Class for a spherical morphology
    SphericalMorphology();
    SphericalMorphology(Morphology *refMorphology, std::shared_ptr<Coordinate3D> position, double radius,
                        std::string description = "");

    void generateSphereRepresentation() final;
    void includeMorphologyXMLOutput(XMLFile *xmlFile, XMLNode *node) final;
    void includeAssociatedMorphologyXMLOutput(XMLFile *xmlFile, XMLNode *node) final;

    std::string generatePovObject();
    double getRadius() { return sphere_rep_.front()->getRadius(); };

    std::string getTypeName() { return "SphericalMorphology"; };

private:
    double radius;
};

#endif    /* SPHERESURFACE_H */

