//  Copyright by Christoph Saffer, Paul Rudolph, Sandra Timme, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "external/doctest/doctest.h"

#include "testUnits.h"
#include <memory>
#include <string>

#include "analyser/pair_measurement.h"
#include "analyser/histogram_measurment.h"
#include "analyser/InSituMeasurements.h"


TEST_CASE ("Check Pair Measurements") {
    const auto measurement = std::make_unique<PairMeasurement>("TEST", "Case1", "Case2");
    measurement->addValuePairs(1, 2);
    std::ostringstream out;
    out << *measurement;
    CHECK(out.str() == "TEST;1;2;\n");
    measurement->cache["Case2"] = 5;
    measurement->addValuesFromCache(1, "Case2");
    out.str("");
    out << *measurement;
    CHECK(out.str() == "TEST;1;2;\nTEST;1;5;\n");
    measurement->cache["Case2"] = std::vector<int>{5, 2};
    measurement->addValuesFromCache(3, "Case2");
    out.str("");
    out << *measurement;
    CHECK(out.str() == "TEST;1;2;\nTEST;1;5;\nTEST;3;5,2,;\n");
}

TEST_CASE ("Check Histogram Measurements") {
    const auto measurement = std::make_unique<HistogramMeasurement>("TEST", "Case1", "Case2");
    measurement->addValues(std::make_pair("Case1", 1), std::make_pair("Case2", 2));
    measurement->addValues(std::make_pair("Case2", 5), std::make_pair("Case2", 2));
    std::ostringstream out;
    out << *measurement;
    CHECK(out.str() == "TEST;1;2;\nTEST;;5;\nTEST;;2;\n");
    out.str("");
    measurement->cache["Case1"] = 2;
    measurement->addValuesFromCache("Case1");
    out << *measurement;
    CHECK(out.str() == "TEST;1;2;\nTEST;2;5;\nTEST;;2;\n");
}

TEST_CASE("Check distance function") {
    Coordinate3D test1{1.0, 2.0, 4.0};
    Coordinate3D test2{4.0, 6.0, 4.0};
    double result = 5.0;
    CHECK(test1.calculateEuclidianDistance(test2) == result);
}

// SphericCoordinate3D.cpp
TEST_CASE("Check distance between two spherical coordinates") {
    double r1 = 1.5;
    double r2 = 78.2;
    double distance, result;
    SphericCoordinate3D sphCoord1 = {r1, 0.5 * M_PI, 0.5};
    SphericCoordinate3D sphCoord2 = {r2, 0.5 * M_PI, 0.5};
    SphericCoordinate3D sphCoord3 = {r2, 0.5 * M_PI, 0.5 + M_PI};

    SphericCoordinate3D sphCoord4 = {r2, 0, 0.3};
    SphericCoordinate3D sphCoord5 = {r2, M_PI, 5.31};

    distance = sphCoord2.calculateEuclidianDistance(sphCoord1);
    result = r2 - r1;
    CHECK(distance == result);

    distance = sphCoord3.calculateEuclidianDistance(sphCoord2);
    result = 2 * r2;
    CHECK(distance == result);

    distance = sphCoord5.calculateEuclidianDistance(sphCoord4);
    result = 2 * r2;
    CHECK(distance == result);
}

TEST_CASE("Check transformation between cartesian and spherical coordinates") {
    SphericCoordinate3D sphCoord0 = {23.68, 1.7, 4.33};
    SphericCoordinate3D sphCoord1 = {23.68, 1.7, 4.33};
    Coordinate3D carCoord1;
    int n = 5; // transform n times between spherical and cartesian coordinates
    int r_increase = 1.5;

    for (int i = 0; i < n; ++i) {
        carCoord1 = abm::util::toCartesianCoordinates(sphCoord1);
        sphCoord1 = abm::util::toSphericCoordinates(carCoord1);
        sphCoord1.r += r_increase;
    }

    double distance = sphCoord0.calculateEuclidianDistance(sphCoord1);
    double result = n * r_increase;
    CHECK(doctest::Approx(distance) == result);
}

TEST_CASE("Check distance on spheres") {
    double r = 95.0;
    SphericCoordinate3D sphCoord1 = {r, 0.5 * M_PI, 0.3};
    SphericCoordinate3D sphCoord2 = {r, 0.5 * M_PI, 0.3 + 0.5 * M_PI};

    double distance = sphCoord1.calculateSphericalDistance(sphCoord2);
    double result = 2 * r * M_PI * 0.25;
    CHECK(distance == result);

    SphericCoordinate3D sphCoord3 = {r, M_PI, 0.3};
    SphericCoordinate3D sphCoord4 = {r, 0, 4.43};

    distance = sphCoord3.calculateSphericalDistance(sphCoord4);
    result = 2 * r * M_PI * 0.5;
    CHECK(distance == result);
}
