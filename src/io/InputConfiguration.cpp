//  Copyright by Christoph Saffer, Paul Rudolph, Sandra Timme, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#include <string>

#include "io/InputConfiguration.h"

InputConfiguration::InputConfiguration(const std::string &configFilename) {
    xRootNode = XMLNode::openFileHelper(configFilename.c_str(), "Agent-Based-Framework");
    scanXMLGraphForNodeIds();
}

InputConfiguration::InputConfiguration(const std::string &configFilename, const std::string &rootTag) {
    xRootNode = XMLNode::openFileHelper(configFilename.c_str(), rootTag.c_str());
    scanXMLGraphForNodeIds();
}

XMLNode InputConfiguration::getRootNode() const {
    return xRootNode;
}

void InputConfiguration::scanXMLGraphForNodeIds() {
    scanChildNodesRecursiveForIds(xRootNode);
}

void InputConfiguration::scanChildNodesRecursiveForIds(const XMLNode &node) {
    int attributes = node.nAttribute();
    for (int i = 0; i < attributes; i++) {
        XMLAttribute curAttr = node.getAttribute(i);
        std::string name = curAttr.lpszName;
        std::string value = curAttr.lpszValue;
        if (name == "id") {
            nodesById[value] = node;
        }
    }

    int childNodes = node.nChildNode();
    for (int i = 0; i < childNodes; i++) {
        scanChildNodesRecursiveForIds(node.getChildNode(i));
    }
}

const char *InputConfiguration::getValueOfXMLNodeId(const std::string &id) {
    const char *retVal = "";
    XMLNode valueNode = getXMLNodebyId(id);
    std::string nameOfNode = valueNode.getName();
    if (!valueNode.isEmpty()) {
        if (nameOfNode == "Value") {
            retVal = valueNode.getAttribute("value");
        }
    }

    return retVal;
}

std::string InputConfiguration::getStringDataFieldByName(XMLNode *currentNode, const std::string &name) {
    std::string returnVal;
    int nodes = currentNode->nChildNode("DataField");
    for (int i = 0; i < nodes; i++) {
        std::string currentString = currentNode->getChildNode("DataField", i).getAttribute("name");
        if (currentString == name) {
            returnVal = currentNode->getChildNode("DataField", i).getChildNode("Value").getAttribute("value");
            break;
        }
    }
    return returnVal;
}

double InputConfiguration::getDoubleDataFieldValueByName(XMLNode *currentNode, const std::string &name) {
    double returnVal = 0;
    int nodes = currentNode->nChildNode("DataField");
    for (int i = 0; i < nodes; i++) {
        std::string currentString = currentNode->getChildNode("DataField", i).getAttribute("name");
        if (currentString == name) {
            returnVal = atof(currentNode->getChildNode("DataField", i).getChildNode("Value").getAttribute("value"));
            break;
        }
    }
    return returnVal;
}

double InputConfiguration::getDoubleDataFieldStddevValueByName(XMLNode *currentNode, const std::string &name) {
    double returnVal = 0;
    int nodes = currentNode->nChildNode("DataField");
    for (int i = 0; i < nodes; i++) {
        std::string currentString = currentNode->getChildNode("DataField", i).getAttribute("name");
        if (currentString == name) {
            returnVal = atof(currentNode->getChildNode("DataField", i).getChildNode("Value").getAttribute("stddev"));
            break;
        }
    }
    return returnVal;
}

void InputConfiguration::setDoubleDataFieldValueByName(XMLNode *currentNode, const std::string &name, double newValue) {
    int nodes = currentNode->nChildNode("DataField");
    for (int i = 0; i < nodes; i++) {
        std::string currentString = currentNode->getChildNode("DataField", i).getAttribute("name");
        if (currentString == name) {
            XMLNode valueNode = currentNode->getChildNode("DataField", i).getChildNode("Value");
            valueNode.deleteAttribute("value");
            std::ostringstream dval;
            dval << newValue;
            valueNode.addAttribute("value", dval.str().c_str());
            break;
        }
    }
}

unsigned int InputConfiguration::getUIntDataFieldValueByName(XMLNode *currentNode, const std::string &name) {
    unsigned int returnVal = 0;
    int nodes = currentNode->nChildNode("DataField");
    for (int i = 0; i < nodes; i++) {
        std::string currentString = currentNode->getChildNode("DataField", i).getAttribute("name");
        if (currentString == name) {
            returnVal = (unsigned int) atoi(
                    currentNode->getChildNode("DataField", i).getChildNode("Value").getAttribute("value"));
            break;
        }
    }
    return returnVal;
}

int InputConfiguration::getIntDataFieldValueByName(XMLNode *currentNode, const std::string &name) {
    int returnVal = 0;
    int nodes = currentNode->nChildNode("DataField");
    for (int i = 0; i < nodes; i++) {
        std::string currentString = currentNode->getChildNode("DataField", i).getAttribute("name");
        if (currentString == name) {
            returnVal = (int) atoi(
                    currentNode->getChildNode("DataField", i).getChildNode("Value").getAttribute("value"));
            break;
        }
    }
    return returnVal;
}

Coordinate3D InputConfiguration::getCoordinateDataFieldValueByName(XMLNode *currentNode, const std::string &name) {
    double x = 0, y = 0, z = 0;
    int nodes = currentNode->nChildNode("DataField");
    for (int i = 0; i < nodes; i++) {
        std::string currentString = currentNode->getChildNode("DataField", i).getAttribute("name");
        //if(currentNode->getChildNode("DataField", i).)
        if (currentString == name) {
            XMLNode coordinateValueNode = currentNode->getChildNode("DataField", i).getChildNode("Values");
            if (!coordinateValueNode.isEmpty()) {
                x = atof(coordinateValueNode.getAttribute("xValue"));
                y = atof(coordinateValueNode.getAttribute("yValue"));
                z = atof(coordinateValueNode.getAttribute("zValue"));
            } else {
                coordinateValueNode = currentNode->getChildNode("DataField", i).getChildNode("Value");
                x = atof(coordinateValueNode.getAttribute("xValue"));
                y = atof(coordinateValueNode.getAttribute("yValue"));
                z = atof(coordinateValueNode.getAttribute("zValue"));
            }
            break;
        }
    }
    return Coordinate3D{x, y, z};
}

XMLNode InputConfiguration::getValueNodeByDataFieldName(XMLNode *currentNode, const std::string &name) {
    XMLNode returnNode;
    int nodes = currentNode->nChildNode("DataField");
    for (int i = 0; i < nodes; i++) {
        std::string currentString = currentNode->getChildNode("DataField", i).getAttribute("name");
        if (currentString == name) {
            returnNode = currentNode->getChildNode("DataField", i).getChildNode("Values");
            break;
        }
    }
    return returnNode;
}

void InputConfiguration::writeFile(const std::string &filepathname) const {
    xRootNode.writeToFile(filepathname.c_str(), "utf-8", '1');
}
SphericCoordinate3D InputConfiguration::getSphericVectorDataFieldByName(XMLNode *currentNode, const std::string &name) {
    double rVal = 0, phiVal = 0, thetaVal = 0;
    int nodes = currentNode->nChildNode("DataField");
    for (int i = 0; i < nodes; i++) {
        std::string currentString = currentNode->getChildNode("DataField", i).getAttribute("name");
        if (currentString == name) {
            XMLNode valueNode = currentNode->getChildNode("DataField", i).getChildNode("Value");
            rVal = atof(valueNode.getAttribute("rValue"));
            phiVal = atof(valueNode.getAttribute("phiValue"));
            thetaVal = atof(valueNode.getAttribute("thetaValue"));
            break;
        }
    }
    return SphericCoordinate3D{rVal, thetaVal, phiVal};
}
ColorRGB InputConfiguration::getColorDataFieldByName(XMLNode *currentNode, const std::string &name) {
    double r = 0, g = 0, b = 0, t = 0;
    int nodes = currentNode->nChildNode("DataField");
    for (int i = 0; i < nodes; i++) {
        std::string currentString = currentNode->getChildNode("DataField", i).getAttribute("name");
        if (currentString == name) {
            XMLNode valueNode = currentNode->getChildNode("DataField", i).getChildNode("Value");
            r = atof(valueNode.getAttribute("red"));
            g = atof(valueNode.getAttribute("green"));
            b = atof(valueNode.getAttribute("blue"));
            t = atof(valueNode.getAttribute("transmit"));
            break;
        }
    }
    return ColorRGB(r, g, b, t);
}
