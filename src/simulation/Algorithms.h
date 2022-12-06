//  Copyright by Christoph Saffer, Paul Rudolph, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#ifndef ALGORITHMS_H
#define    ALGORITHMS_H

#include <boost/numeric/ublas/matrix.hpp>

#include "basic/Randomizer.h"
#include "basic/SphericCoordinate3D.h"

class Algorithms {
public:
  // Class for providing useful utility functions.
    static std::vector<unsigned int> generateRandomPermutation(Randomizer *randomizer, unsigned int size);

    static double nChoosek(unsigned long n, unsigned long k);
    static double bernoulliProbability(unsigned long n, unsigned long k, double p);
    static double retrieveDirectionAngleAlpha(SphericCoordinate3D ownPos, SphericCoordinate3D goalPos);

    static double interpolateBilinearOnTriangle(Coordinate3D pointOfInterest, Coordinate3D c1, double w1, Coordinate3D c2, double w2,
                                  Coordinate3D c3, double w3);
    static double interpolateBiquadraticOnTriangle(Coordinate3D pointOfInterest, Coordinate3D c1, double w1, Coordinate3D gradC1,
                                     Coordinate3D c2, double w2, Coordinate3D gradC2, Coordinate3D c3, double w3,
                                     Coordinate3D gradC3);
    static double meanAngle(double *angles, int size);
    bool InvertMatrix(const boost::numeric::ublas::matrix<double> &input, boost::numeric::ublas::matrix<double> &inverse);

private:
    static void swap(std::vector<unsigned int> *vectorToSwap, unsigned int i, unsigned int j);

};

#endif    /* ALGORITHMS_H */

