//  Copyright by Christoph Saffer, Paul Rudolph, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#include <basic/SphericCoordinate3D.h>
#include <sstream>
#include <vector>

#include "basic/Coordinate3D.h"
#include "visualisation/PovRayObject.h"

std::string PovRayObject::getCamera(Coordinate3D location, Coordinate3D look_at, double angle) {
    std::ostringstream ss;
    ss << "camera {" << '\n';
    ss << "\tlocation <" << location.x << "," << location.y << "," << location.z << ">\n";
    ss << "\tlook_at <" << look_at.x << "," << look_at.y << "," << look_at.z << ">\n";
    if (angle != 53.0) {
        ss << "\tangle " << angle << '\n';
    }
    ss << "\tright x*image_width/image_height" << '\n';
    ss << "\t }" << '\n';

    return ss.str();

}

std::string PovRayObject::getBackground(const ColorRGB &color) {
    std::ostringstream ss;
    ss << "background {" << '\n';
    ss << "\tcolor rgb <" << color.getRed() << ",\t" << color.getGreen() << ",\t" << color.getBlue() << ">\t}" << '\n';
    return ss.str();

}

std::string
PovRayObject::getCylinder(Coordinate3D firstCenter, Coordinate3D secondCenter, double radius, ColorRGB color) {
    std::ostringstream ss;
    ss << "cylinder {" << '\n';
    ss << "\t<" << firstCenter.x << "," << firstCenter.y << "," << firstCenter.z << ">,<" << secondCenter.x << ","
       << secondCenter.y << "," << secondCenter.z << ">," << radius << '\n';
    ss << "\tpigment { color " << color.printPovColorRGBT() << " }" << '\n';
    ss << "\t}" << '\n';
    return ss.str();
}

std::string PovRayObject::getLightsource(Coordinate3D position) {
    std::ostringstream ss;
    ss << "light_source {" << '\n';
    ss << "\t<" << position.x << "," << position.y << "," << position.z << ">\n";
    ss << "\trgb <1,1,1>" << '\n';
    ss << "\tshadowless" << '\n';
    ss << "\t}" << '\n';

    return ss.str();
}

std::string PovRayObject::getSphere(Coordinate3D center, double radius, ColorRGB color) {
    std::ostringstream ss;
    ss << "sphere {" << '\n';
    ss << "\t<" << center.x << "," << center.y << "," << center.z << ">," << "\n";
    ss << "\t" << radius << '\n';
    ss << "\t" << "pigment { color " << color.printPovColorRGBT() << " }" << '\n';
    ss << "\t}" << '\n';

    return ss.str();
}

std::string PovRayObject::getBox(Coordinate3D leftDown, Coordinate3D rightUp, ColorRGB color) {
    std::ostringstream ss;
    ss << "box {" << '\n';
    ss << "\t<" << leftDown.x << "," << leftDown.y << "," << leftDown.z << ">,<" << rightUp.x << "," << rightUp.y << ","
       << rightUp.z << ">\n";
    ss << "\t" << "pigment { color " << color.printPovColorRGBT() << " }" << '\n';
    ss << "\t}" << '\n';

    return ss.str();
}

std::string
PovRayObject::getBox(Coordinate3D left_down, Coordinate3D right_up, ColorRGB color, SphericCoordinate3D rotate,
                     Coordinate3D translate) {
    std::ostringstream ss;
    ss << "box {" << '\n';
    ss << "\t<" << left_down.x << "," << left_down.y << "," << left_down.z << ">,<" << right_up.x << "," << right_up.y
       << "," << right_up.z << ">\n";
    ss << "\t rotate <0," << rotate.theta * 180.0 / M_PI << ",0>" << '\n';
    ss << "\t rotate <0,0," << rotate.phi * 180.0 / M_PI << ">" << '\n';
    ss << "\t translate <" << translate.x << "," << translate.y << "," << translate.z << ">\n";
    ss << "\t" << "pigment { color " << color.printPovColorRGBT() << " }" << '\n';
    ss << "\t}" << '\n';

    return ss.str();
}

std::string
PovRayObject::getTorus(Coordinate3D center, double innerRad, double outerRad, Coordinate3D rotation, ColorRGB *color) {
    std::ostringstream oss;
    oss << "torus {\n";
    oss << "\t" << outerRad << ", " << innerRad << "\n";
    //important: rotation before translation
    oss << "\trotate <" << rotation.x << "," << rotation.y << "," << rotation.z << ">\n";
    oss << "\t translate <" << center.x << "," << center.y << "," << center.z << ">\n";
    oss << "\t" << "pigment { color " << color->printPovColorRGBT() << " }" << '\n';
    oss << "\n}\n";
    return oss.str();
}

std::string PovRayObject::getIntersection(const std::string &objOne, const std::string &objTwo) {
    std::ostringstream ss;
    ss << "intersection {" << '\n';
    ss << "\t" << objOne << '\n';
    ss << "\t" << objTwo << '\n';
    ss << "\t}" << '\n';

    return ss.str();
}

std::string PovRayObject::getIntersection(const std::string &objOne, const std::vector<std::string> &objectsTwo) {
    std::ostringstream ss;
    ss << "intersection {" << '\n';
    ss << "\t" << objOne << '\n';

    for (const auto &it : objectsTwo) {
        ss << "\t" << it << '\n';
    }

    ss << "\t}" << '\n';

    return ss.str();
}

std::string PovRayObject::getDifference(const std::string &objOne, const std::string &objTwo) {
    std::ostringstream ss;
    ss << "difference {" << '\n';
    ss << "\t" << objOne << '\n';
    ss << "\t" << objTwo << '\n';
    ss << "\t}" << '\n';

    return ss.str();
}

std::string PovRayObject::getDifference(const std::string &objOne, const std::vector<std::string> &objectsTwo) {
    std::ostringstream ss;
    ss << "difference {" << '\n';
    ss << "\t" << objOne << '\n';
    for (const auto &it : objectsTwo) {
        ss << "\t" << it << '\n';
    }
    ss << "\t}" << '\n';

    return ss.str();
}