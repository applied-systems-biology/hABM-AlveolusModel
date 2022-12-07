//  Copyright by Christoph Saffer, Paul Rudolph, Sandra Timme, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#ifndef SPHERICCOORDINATES3D_H
#define    SPHERICCOORDINATES3D_H

#include <sstream>

struct SphericCoordinate3D {
    double r;
    double theta;
    double phi;

    SphericCoordinate3D operator+(const SphericCoordinate3D &vec) const;
    SphericCoordinate3D &operator+=(const SphericCoordinate3D &vec);
    SphericCoordinate3D operator-(const SphericCoordinate3D &vec) const;
    SphericCoordinate3D &operator-=(const SphericCoordinate3D &vec);
    SphericCoordinate3D operator*(double value) const;
    SphericCoordinate3D &operator*=(double value);

    void setAntipode();

    [[nodiscard]] std::string printCoordinates() const;
    [[nodiscard]] double calculateSphericalDistance(const SphericCoordinate3D &s_cood3d) const noexcept;
    [[nodiscard]] double calculateEuclidianDistance(const SphericCoordinate3D &s_cood3d) const noexcept;
};

#endif    /* SPHERICCOORDINATES3D_H */

