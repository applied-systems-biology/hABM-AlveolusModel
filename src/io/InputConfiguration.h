//  Copyright by Christoph Saffer, Paul Rudolph, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#ifndef INPUTCONFIGURATION_H
#define    INPUTCONFIGURATION_H

#include <sstream>
#include <map>

#include "external/xmlParser/xmlParser.h"
#include "basic/Coordinate3D.h"
#include "basic/SphericCoordinate3D.h"
#include "basic/ColorRGB.h"

class InputConfiguration {
public:
  // Class for reading input files
    InputConfiguration() = default;
    explicit InputConfiguration(const std::string &configFilename);
    InputConfiguration(const std::string &configFilename, const std::string &rootTag);

    [[nodiscard]] XMLNode getRootNode() const;
    void scanXMLGraphForNodeIds();
    void scanChildNodesRecursiveForIds(const XMLNode &node);
    [[nodiscard]] XMLNode getXMLNodebyId(const std::string &id) const {
        if (nodesById.find(id) != nodesById.end())
            return nodesById.at(id);
        else return XMLNode();
    };
    const char *getValueOfXMLNodeId(const std::string &id);
    static void setDoubleDataFieldValueByName(XMLNode *, const std::string &name, double newValue);
    void writeFile(const std::string &filepathname) const;

    static std::string getStringDataFieldByName(XMLNode *, const std::string &name);
    static double getDoubleDataFieldValueByName(XMLNode *, const std::string &name);
    static double getDoubleDataFieldStddevValueByName(XMLNode *, const std::string &name);
    static Coordinate3D getCoordinateDataFieldValueByName(XMLNode *, const std::string &name);
    static unsigned int getUIntDataFieldValueByName(XMLNode *, const std::string &name);
    static int getIntDataFieldValueByName(XMLNode *, const std::string &name);
    static XMLNode getValueNodeByDataFieldName(XMLNode *, const std::string &name);
    static SphericCoordinate3D getSphericVectorDataFieldByName(XMLNode *currentNode, const std::string &name);
    static ColorRGB getColorDataFieldByName(XMLNode *currentNode, const std::string &name);

private:
    XMLNode xRootNode;
    std::map<std::string, XMLNode> nodesById;
};

#endif    /* INPUTCONFIGURATION_H */

