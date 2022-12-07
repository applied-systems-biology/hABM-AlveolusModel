//  Copyright by Christoph Saffer, Paul Rudolph, Sandra Timme, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#ifndef POVRAYOBJECT_H
#define    POVRAYOBJECT_H

#include "basic/Coordinate3D.h"
#include "basic/ColorRGB.h"

class PovRayObject {
public:
    static std::string getSphere(Coordinate3D center, double radius, ColorRGB color);
    static std::string getBox(Coordinate3D leftDown, Coordinate3D rightUp, ColorRGB color);
    static std::string getBox(Coordinate3D leftDown, Coordinate3D rightUp, ColorRGB color, SphericCoordinate3D rotate,
                              Coordinate3D translate);
    static std::string getCylinder(Coordinate3D firstCenter, Coordinate3D secondCenter, double radius, ColorRGB color);
    static std::string getLightsource(Coordinate3D);
    static std::string getBackground(const ColorRGB &color);
    static std::string getCamera(Coordinate3D location, Coordinate3D look_at, double angle = 53.0);
    static std::string
    getTorus(Coordinate3D center, double innerRad, double outerRad, Coordinate3D rotation, ColorRGB *color);
    static std::string getIntersection(const std::string &objOne, const std::string &objTwostatic);
    static std::string getIntersection(const std::string &objOne, const std::vector<std::string> &objectsTwo);
    static std::string getDifference(const std::string &objOne, const std::string &objTwo);
    static std::string getDifference(const std::string &objOne, const std::vector<std::string> &objectsTwo);
};

#endif    /* POVRAYOBJECT_H */

