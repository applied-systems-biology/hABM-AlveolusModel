//  Copyright by Christoph Saffer, Paul Rudolph, Sandra Timme, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#include <cmath>
#include <utils/misc_util.h>

#include "basic/Randomizer.h"
#include "basic/SphericCoordinate3D.h"
#include "utils/macros.h"

double SphericCoordinate3D::calculateSphericalDistance(const SphericCoordinate3D &s_cood3d) const noexcept {
    if (std::abs(r - s_cood3d.r) > r / 10.0) {
        ERROR_STDERR("Warning: Spherical coordinates are not on the same surface (differ in radius).");
        exit(1);
    }
    double dphi = s_cood3d.phi - phi;
    double dalpha = acos(sin(theta) * sin(s_cood3d.theta) * cos(dphi) + cos(theta) * cos(s_cood3d.theta));
    return dalpha * r;
}

double SphericCoordinate3D::calculateEuclidianDistance(const SphericCoordinate3D &s_cood3d) const noexcept {
    return abm::util::toCartesianCoordinates(*this).calculateEuclidianDistance(
            abm::util::toCartesianCoordinates(s_cood3d));
}

void SphericCoordinate3D::setAntipode() {
    theta = M_PI - theta;
    if (phi >= M_PI) {
        phi -= M_PI;
    } else {
        phi += M_PI;
    }
}

std::string SphericCoordinate3D::printCoordinates() const {
    std::stringstream ss;
    ss << "(" << r << ", " << theta << ", " << phi << ")";
    return ss.str();
}

SphericCoordinate3D SphericCoordinate3D::operator+(const SphericCoordinate3D &vec) const {
    return abm::util::toSphericCoordinates(
            abm::util::toCartesianCoordinates(*this) + abm::util::toCartesianCoordinates(vec));
}

SphericCoordinate3D &SphericCoordinate3D::operator+=(const SphericCoordinate3D &vec) {
    const auto tmp_value = abm::util::toSphericCoordinates(
            abm::util::toCartesianCoordinates(*this) + abm::util::toCartesianCoordinates(vec));
    r = tmp_value.r;
    theta = tmp_value.theta;
    phi = tmp_value.phi;
    return *this;
}

SphericCoordinate3D &SphericCoordinate3D::operator-=(const SphericCoordinate3D &vec) {
    const auto tmp_value = abm::util::toSphericCoordinates(
            abm::util::toCartesianCoordinates(*this) - abm::util::toCartesianCoordinates(vec));
    r = tmp_value.r;
    theta = tmp_value.theta;
    phi = tmp_value.phi;
    return *this;
}

SphericCoordinate3D SphericCoordinate3D::operator-(const SphericCoordinate3D &vec) const {
    return abm::util::toSphericCoordinates(
            abm::util::toCartesianCoordinates(*this) - abm::util::toCartesianCoordinates(vec));
}

SphericCoordinate3D SphericCoordinate3D::operator*(double value) const {
    return abm::util::toSphericCoordinates(abm::util::toCartesianCoordinates(*this) * value);
}

SphericCoordinate3D &SphericCoordinate3D::operator*=(double value) {
    const auto tmp_value = abm::util::toSphericCoordinates(abm::util::toCartesianCoordinates(*this) * value);
    r = tmp_value.r;
    theta = tmp_value.theta;
    phi = tmp_value.phi;
    return *this;
}
