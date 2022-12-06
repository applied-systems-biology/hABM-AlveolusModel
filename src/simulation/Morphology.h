//  Copyright by Christoph Saffer, Paul Rudolph, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#ifndef SURFACE_H
#define    SURFACE_H

#include <list>

#include "basic/ColorRGB.h"
#include "basic/Coordinate3D.h"
#include "io/XMLFile.h"
#include "simulation/morphology/MorphologyElement.h"

class Cell;

class Morphology {
public:
  // Class for handling the morphology of cells
    Morphology() = default;

    Morphology(const std::string &color, Cell *cell);
    virtual ~Morphology() = default;
    virtual std::string getTypeName();
    std::string generatePovObject();
    double getMaximumDistanceFromCenter();
    void includeMorphologyXMLOutput(XMLFile *xmlFile, XMLNode *node);
    ColorRGB *getColorRGB();
    void setColorRGB(std::unique_ptr<ColorRGB> col);
    Cell *getCellThisBelongsTo() { return cell_this_belongs_to_; };
    void appendAssociatedCellpart(std::unique_ptr<MorphologyElement> morphElement);
    std::vector<SphereRepresentation *> getAllSpheresOfThis();
    SphereRepresentation *getBasicSphereOfThis();
    double getVolume();

protected:
    std::unique_ptr<ColorRGB> color_rgb_{};
    std::list<std::unique_ptr<MorphologyElement>> morphologyElements;
    Cell *cell_this_belongs_to_{};

};

#endif    /* SURFACE_H */

