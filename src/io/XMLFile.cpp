//  Copyright by Christoph Saffer, Paul Rudolph, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#include <iostream>
#include <sstream>
#include <string>

#include "io/XMLFile.h"
#include "utils/macros.h"

XMLFile::XMLFile() {
    rootNode = XMLNode::createXMLTopNode("Agent-Based-Framework");
    filename = "";
    includeHeader();
}

XMLFile::XMLFile(std::string nameOfRootNode) {
    rootNode = XMLNode::createXMLTopNode(nameOfRootNode.c_str());
    filename = "";
    includeHeader();
}

XMLFile::XMLFile(std::string filename, std::string nameOfRootNode) {
    rootNode = XMLNode::openFileHelper(filename.c_str(), nameOfRootNode.c_str());
    this->filename = filename;
}

void XMLFile::includeHeader() {
    XMLNode header = addChildToRootNode("Header");
    header.addAttribute("copyright", "Johannes Pollmächer");
}

XMLNode XMLFile::getRootNode() {
    return rootNode;
}

XMLNode XMLFile::addChildToRootNode(std::string tag) {
    return rootNode.addChild(tag.c_str());
}

XMLNode XMLFile::addChildToNode(XMLNode xmlNode, std::string tag) {
    return xmlNode.addChild(tag.c_str());
}

void XMLFile::addDataFieldToNode(XMLNode xmlNode, std::string datename, std::string conversion, std::string dataType,
                                 std::string value) {
    XMLNode dF = xmlNode.addChild("DataField");
    dF.addAttribute("name", datename.c_str());
    dF.addAttribute("conversion", conversion.c_str());
    dF.addAttribute("dataType", dataType.c_str());

    XMLNode dFVal = dF.addChild("Value");
    dFVal.addAttribute("value", value.c_str());
}

void XMLFile::addDataFieldToNode(XMLNode xmlNode, const Coordinate3D &position) {
    XMLNode dF = xmlNode.addChild("DataField");
    dF.addAttribute("name", "position");
    dF.addAttribute("conversion", "continuous");
    dF.addAttribute("dataType", "double");

    XMLNode dFVal = dF.addChild("Value");
    std::ostringstream ssX;
    std::ostringstream ssY;
    std::ostringstream ssZ;
    ssX << position.x;
    dFVal.addAttribute("xValue", ssX.str().c_str());
    ssY << position.y;
    dFVal.addAttribute("yValue", ssY.str().c_str());
    ssZ << position.z;
    dFVal.addAttribute("zValue", ssZ.str().c_str());
}

void XMLFile::addDataFieldToNode(XMLNode xmlNode, const Coordinate3D &position, const std::string &fieldName) {
    XMLNode dF = xmlNode.addChild("DataField");
    dF.addAttribute("name", fieldName.c_str());
    dF.addAttribute("conversion", "continuous");
    dF.addAttribute("dataType", "double");

    XMLNode dFVal = dF.addChild("Value");
    std::ostringstream ssX, ssY, ssZ;
    ssX << position.x;
    dFVal.addAttribute("xValue", ssX.str().c_str());
    ssY << position.y;
    dFVal.addAttribute("yValue", ssY.str().c_str());
    ssZ << position.z;
    dFVal.addAttribute("zValue", ssZ.str().c_str());
}

void XMLFile::addDataFieldToNode(XMLNode xmlNode, const SphericCoordinate3D &position) {
    XMLNode dF = xmlNode.addChild("DataField");
    dF.addAttribute("name", "position");
    dF.addAttribute("conversion", "continuous");
    dF.addAttribute("dataType", "double");

    XMLNode dFVal = dF.addChild("Value");
    std::ostringstream ssX, ssY, ssZ;
    ssX << position.r;
    dFVal.addAttribute("rValue", ssX.str().c_str());
    ssY << position.phi;
    dFVal.addAttribute("phiValue", ssY.str().c_str());
    ssZ << position.theta;
    dFVal.addAttribute("thetaValue", ssZ.str().c_str());
}

void XMLFile::addDataFieldToNode(XMLNode xmlNode, const SphericCoordinate3D &position, const std::string &fieldName) {
    XMLNode dF = xmlNode.addChild("DataField");
    dF.addAttribute("name", fieldName.c_str());
    dF.addAttribute("conversion", "continuous");
    dF.addAttribute("dataType", "double");

    XMLNode dFVal = dF.addChild("Value");
    std::ostringstream ssX, ssY, ssZ;
    ssX << position.r;
    dFVal.addAttribute("rValue", ssX.str().c_str());
    ssY << position.phi;
    dFVal.addAttribute("phiValue", ssY.str().c_str());
    ssZ << position.theta;
    dFVal.addAttribute("thetaValue", ssZ.str().c_str());
}

void XMLFile::addDataFieldToNode(XMLNode xmlNode, const ColorRGB &colorRGB) {
    XMLNode dF = xmlNode.addChild("DataField");
    dF.addAttribute("name", "color");
    dF.addAttribute("conversion", "continuous");
    dF.addAttribute("dataType", "rgbt");

    XMLNode dFVal = dF.addChild("Value");
    std::ostringstream ssR, ssG, ssB, ssT;
    ssR << colorRGB.getRed();
    dFVal.addAttribute("red", ssR.str().c_str());
    ssG << colorRGB.getGreen();
    dFVal.addAttribute("green", ssG.str().c_str());
    ssB << colorRGB.getBlue();
    dFVal.addAttribute("blue", ssB.str().c_str());
    ssT << colorRGB.getTransmit();
    dFVal.addAttribute("transmit", ssT.str().c_str());
}

void XMLFile::setFilename(std::string filename) {
    this->filename = filename;
}

void XMLFile::writeFile() {
    if (filename.length() > 0) {
        rootNode.writeToFile(filename.c_str(), "utf-8", '1');
    } else {
        ERROR_STDERR("no output possible... no file specified");
    }
}

void XMLFile::writeFile(std::string filepathname) {
    rootNode.writeToFile(filepathname.c_str(), "utf-8", '1');
}

Coordinate3D XMLFile::getCoordinateDataFieldByName(XMLNode *currentNode, std::string name) {
    double xVal = 0, yVal = 0, zVal = 0;
    int nodes = currentNode->nChildNode("DataField");
    for (int i = 0; i < nodes; i++) {
        std::string currentString = currentNode->getChildNode("DataField", i).getAttribute("name");
        if (currentString.compare(name) == 0) {
            XMLNode valueNode = currentNode->getChildNode("DataField", i).getChildNode("Value");
            xVal = atof(valueNode.getAttribute("xValue"));
            yVal = atof(valueNode.getAttribute("yValue"));
            zVal = atof(valueNode.getAttribute("zValue"));
            break;
        }
    }
    return Coordinate3D{xVal, yVal, zVal};
}

SphericCoordinate3D XMLFile::getSphericVectorDataFieldByName(XMLNode *currentNode, std::string name) {
    double rVal = 0, phiVal = 0, thetaVal = 0;
    int nodes = currentNode->nChildNode("DataField");
    for (int i = 0; i < nodes; i++) {
        std::string currentString = currentNode->getChildNode("DataField", i).getAttribute("name");
        if (currentString.compare(name) == 0) {
            XMLNode valueNode = currentNode->getChildNode("DataField", i).getChildNode("Value");
            rVal = atof(valueNode.getAttribute("rValue"));
            phiVal = atof(valueNode.getAttribute("phiValue"));
            thetaVal = atof(valueNode.getAttribute("thetaValue"));
            break;
        }
    }
    return SphericCoordinate3D{rVal, thetaVal, phiVal};
}

ColorRGB XMLFile::getColorDataFieldByName(XMLNode *currentNode, std::string name) {
    double r = 0, g = 0, b = 0, t = 0;
    int nodes = currentNode->nChildNode("DataField");
    for (int i = 0; i < nodes; i++) {
        std::string currentString = currentNode->getChildNode("DataField", i).getAttribute("name");
        if (currentString.compare(name) == 0) {
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

std::string XMLFile::getStringDataFieldByName(XMLNode *currentNode, std::string name) {
    std::string returnVal = "";
    int nodes = currentNode->nChildNode("DataField");
    for (int i = 0; i < nodes; i++) {
        std::string currentString = currentNode->getChildNode("DataField", i).getAttribute("name");
        if (currentString.compare(name) == 0) {
            returnVal = currentNode->getChildNode("DataField", i).getChildNode("Value").getAttribute("value");
            break;
        }
    }
    return returnVal;
}

double XMLFile::getDoubleDataFieldValueByName(XMLNode *currentNode, std::string name) {
    double returnVal = 0;
    int nodes = currentNode->nChildNode("DataField");
    for (int i = 0; i < nodes; i++) {
        std::string currentString = currentNode->getChildNode("DataField", i).getAttribute("name");
        if (currentString.compare(name) == 0) {
            returnVal = atof(currentNode->getChildNode("DataField", i).getChildNode("Value").getAttribute("value"));
            break;
        }
    }
    return returnVal;
}

void XMLFile::setDoubleDataFieldValueByName(XMLNode *currentNode, std::string name, double newValue) {
    int nodes = currentNode->nChildNode("DataField");
    for (int i = 0; i < nodes; i++) {
        std::string currentString = currentNode->getChildNode("DataField", i).getAttribute("name");
        if (currentString.compare(name) == 0) {
            XMLNode valueNode = currentNode->getChildNode("DataField", i).getChildNode("Value");
            valueNode.deleteAttribute("value");
            std::ostringstream dval;
            dval << newValue;
            valueNode.addAttribute("value", dval.str().c_str());
            break;
        }
    }
}

unsigned int XMLFile::getUIntDataFieldValueByName(XMLNode *currentNode, std::string name) {
    unsigned int returnVal = 0;
    int nodes = currentNode->nChildNode("DataField");
    for (int i = 0; i < nodes; i++) {
        std::string currentString = currentNode->getChildNode("DataField", i).getAttribute("name");
        if (currentString.compare(name) == 0) {
            returnVal = (unsigned int) atoi(
                    currentNode->getChildNode("DataField", i).getChildNode("Value").getAttribute("value"));
            break;
        }
    }
    return returnVal;
}

XMLNode XMLFile::getValueNodeByDataFieldName(XMLNode *currentNode, std::string name) {
    XMLNode returnNode;
    int nodes = currentNode->nChildNode("DataField");
    for (int i = 0; i < nodes; i++) {
        std::string currentString = currentNode->getChildNode("DataField", i).getAttribute("name");
        if (currentString.compare(name) == 0) {
            returnNode = currentNode->getChildNode("DataField", i).getChildNode("Values");
            break;
        }
    }
    return returnNode;
}