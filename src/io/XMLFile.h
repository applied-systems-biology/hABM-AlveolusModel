//  Copyright by Christoph Saffer, Paul Rudolph, Sandra Timme, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#ifndef XMLFILE_H
#define    XMLFILE_H

#include "basic/ColorRGB.h"
#include "basic/Coordinate3D.h"
#include "external/xmlParser/xmlParser.h"
#include "basic/SphericCoordinate3D.h"

class XMLFile {
public:
  // Class for wrapping XMLFile read and write instructions.
    XMLFile();

    explicit XMLFile(std::string nameOfRootNode);
    XMLFile(std::string filename, std::string nameOfRootNode);
    void includeHeader();
    XMLNode getRootNode();
    XMLNode addChildToRootNode(std::string tag);
    XMLNode addChildToNode(XMLNode xmlNode, std::string tag);
    void addDataFieldToNode(XMLNode xmlNode, std::string datename, std::string optype, std::string dataType,
                            std::string value);
    void addDataFieldToNode(XMLNode xmlNode, const Coordinate3D &position);
    void addDataFieldToNode(XMLNode xmlNode, const Coordinate3D &position, const std::string &fieldName);
    void addDataFieldToNode(XMLNode xmlNode, const SphericCoordinate3D &position);
    void addDataFieldToNode(XMLNode xmlNode, const SphericCoordinate3D &position, const std::string &fieldName);
    void addDataFieldToNode(XMLNode xmlNode, const ColorRGB &colorRGB);
    void setFilename(std::string filename);
    void writeFile();
    void writeFile(std::string filepathname);
    Coordinate3D getCoordinateDataFieldByName(XMLNode *currentNode, std::string name);
    SphericCoordinate3D getSphericVectorDataFieldByName(XMLNode *currentNode, std::string name);
    ColorRGB getColorDataFieldByName(XMLNode *currentNode, std::string name);
    std::string getStringDataFieldByName(XMLNode *currentNode, std::string name);
    double getDoubleDataFieldValueByName(XMLNode *currentNode, std::string name);
    void setDoubleDataFieldValueByName(XMLNode *currentNode, std::string name, double newValue);
    unsigned int getUIntDataFieldValueByName(XMLNode *currentNode, std::string name);
    XMLNode getValueNodeByDataFieldName(XMLNode *currentNode, std::string name);

private:
    XMLNode rootNode;
    std::string filename;

};

#endif    /* XMLFILE_H */

