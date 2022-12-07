//  Copyright by Christoph Saffer, Paul Rudolph, Sandra Timme, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#include "basic/Coordinate3D.h"
#include "basic/Randomizer.h"

double Coordinate3D::calculateEuclidianDistance(const Coordinate3D &coord) const noexcept {
    const double dx = x - coord.x;
    const double dy = y - coord.y;
    const double dz = z - coord.z;
    return sqrt(dx * dx + dy * dy + dz * dz);
}

double Coordinate3D::getMagnitude() const {
    double a;
    if (x == 0 && y == 0 && z == 0) {
        a = 0;
    } else {
        a = sqrt(x * x + y * y + z * z);
    }
    return a;
}

void Coordinate3D::setMagnitude(double length) {
    double magnitude = getMagnitude();
    if (magnitude > 0) {
        *this *= length / magnitude;
    }
}

std::string Coordinate3D::printCoordinates() const {
    std::stringstream ss;
    ss << "(" << x << ", " << y << ", " << z << ")";
    return ss.str();
}

Coordinate3D Coordinate3D::operator+(const Coordinate3D &vec) const {
    return {vec.x + x, vec.y + y, vec.z + z};
}

Coordinate3D &Coordinate3D::operator+=(const Coordinate3D &vec) {
    x += vec.x;
    y += vec.y;
    z += vec.z;
    return *this;
}
Coordinate3D Coordinate3D::operator-(const Coordinate3D &vec) const {
    return {x - vec.x, y - vec.y, z - vec.z};
}
Coordinate3D &Coordinate3D::operator-=(const Coordinate3D &vec) {
    x -= vec.x;
    y -= vec.y;
    z -= vec.z;
    return *this;
}

Coordinate3D Coordinate3D::operator*(double value) const {
    return {x * value, y * value, z * value};
}

Coordinate3D &Coordinate3D::operator*=(double value) {
    x *= value;
    y *= value;
    z *= value;
    return *this;
}

double Coordinate3D::scalarProduct(const Coordinate3D &vec) const {
    return x * vec.x + y * vec.y + z * vec.z;
}

Coordinate3D Coordinate3D::crossProduct(const Coordinate3D &vec) const {
    Coordinate3D new_vec{};
    new_vec.x = y * vec.z - z * vec.y;
    new_vec.y = z * vec.x - x * vec.z;
    new_vec.z = x * vec.y - y * vec.x;
    return new_vec;
}

