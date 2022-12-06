//  Copyright by Christoph Saffer, Paul Rudolph, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#include <vector>
#include <cmath>

#include <boost/numeric/ublas/lu.hpp>

#include <utils/misc_util.h>
#include "simulation/Algorithms.h"
#include "utils/macros.h"

std::vector<unsigned int> Algorithms::generateRandomPermutation(Randomizer *randomizer, unsigned int size) {
    std::vector<unsigned int> permutationVector;
    if (size > 0) {
        //initialize the vector
        for (unsigned int k = 0; k < size; k++) {
            permutationVector.push_back(k);
        }
        for (unsigned int k = 0; k < size - 1; k++) {
            swap(&permutationVector, k, randomizer->generateInt(k, size - 1));
        }
    }
    return permutationVector;
}

void Algorithms::swap(std::vector<unsigned int> *vectorToSwap, unsigned int i, unsigned int j) {
    unsigned int temp;
    temp = (*vectorToSwap)[i];
    (*vectorToSwap)[i] = (*vectorToSwap)[j];
    (*vectorToSwap)[j] = temp;
}

double Algorithms::nChoosek(unsigned long n, unsigned long k) {
    if (k * 2 > n) { k = n - k; }

    if (k > n) return 0;
    if (k == 0) return 1;

    double result = 1.0;
    for (unsigned long i = 1; i <= k; ++i) {
        result *= (n - i + 1);
        result /= i;
    }
    return result;
}

double Algorithms::bernoulliProbability(unsigned long n, unsigned long k, double p) {

    return nChoosek(n, k) * pow(p, k) * pow(1 - p, n - k);
}

double Algorithms::retrieveDirectionAngleAlpha(SphericCoordinate3D ownPos, SphericCoordinate3D goalPos) {
    double alpha = 0;

    double dphi = ownPos.phi - goalPos.phi;
    double dphiAbs = fabs(dphi);

    double c = ownPos.theta;
    double b = goalPos.theta;
    double a = acos(sin(ownPos.theta) * sin(goalPos.theta) * cos(dphi) + cos(ownPos.theta) * cos(goalPos.theta));
//    cout << "[AllveoleSite] distance(1,2)=" << a*radius << '\n';
    if (c == 0) {
        alpha = goalPos.phi;
    } else {
        if (b == 0) {
            alpha = M_PI;
        } else {
            //compute beta by cases
            double beta;
            if (ownPos.phi == goalPos.phi) {
                if (ownPos.theta > goalPos.theta) {
                    beta = 0;
                } else {
                    beta = M_PI;
                }
            } else {
                beta = acos((cos(b) - cos(a) * cos(c)) / (sin(a) * sin(c)));
            }

            //compute alpha by cases
            if ((ownPos.phi <= goalPos.phi && dphiAbs <= M_PI) || (ownPos.phi > goalPos.phi && dphiAbs > M_PI)) {
                alpha = M_PI - beta;
            } else {
                alpha = M_PI + beta;
            }

        }
    }

    return alpha;
}

double
Algorithms::interpolateBilinearOnTriangle(Coordinate3D pointOfInterest, Coordinate3D c1, double w1, Coordinate3D c2,
                                          double w2, Coordinate3D c3, double w3) {

    //first step: rotate all points in the xy-layer
    Coordinate3D triVec1(c2 - c1);
    Coordinate3D triVec2(c3 - c1);

    Coordinate3D normVectorTriangle = triVec1.crossProduct(triVec2);
    Coordinate3D normVectorCoord = Coordinate3D{normVectorTriangle.x, normVectorTriangle.y, normVectorTriangle.z};
    SphericCoordinate3D normVectorSpheric = abm::util::toSphericCoordinates(normVectorCoord);
    double phiNorm = normVectorSpheric.phi;
    double thetaNorm = normVectorSpheric.theta;

    Coordinate3D poiRotPre = abm::util::rotateAboutThetaAndPhi(pointOfInterest, phiNorm, thetaNorm);
    Coordinate3D c1Rot = abm::util::rotateAboutThetaAndPhi(c1, phiNorm, thetaNorm);
    Coordinate3D poiRot{poiRotPre.x, poiRotPre.y, c1Rot.z};
    Coordinate3D c2Rot = abm::util::rotateAboutThetaAndPhi(c2, phiNorm, thetaNorm);
    Coordinate3D c3Rot = abm::util::rotateAboutThetaAndPhi(c3, phiNorm, thetaNorm);

    //step two: find intersection of 2->3 and 1->poi
    double D = (c3Rot.y - c2Rot.y) * (poiRot.x - c1Rot.x) - (poiRot.y - c1Rot.y) * (c3Rot.x - c2Rot.x);

    double xS = ((c3Rot.x - c2Rot.x) * (poiRot.x * c1Rot.y - c1Rot.x * poiRot.y) -
                 (poiRot.x - c1Rot.x) * (c3Rot.x * c2Rot.y - c2Rot.x * c3Rot.y)) / D;
    double yS = ((c1Rot.y - poiRot.y) * (c3Rot.x * c2Rot.y - c2Rot.x * c3Rot.y) -
                 (c2Rot.y - c3Rot.y) * (poiRot.x * c1Rot.y - c1Rot.x * poiRot.y)) / D;
    Coordinate3D rS = Coordinate3D{xS, yS, c1Rot.z};

    //step three: interpolate the concentration bilinearly
    double d23 = c2Rot.calculateEuclidianDistance(c3Rot);
    double d2S = c2Rot.calculateEuclidianDistance(rS);
    double d1S = c1Rot.calculateEuclidianDistance(rS);
    double d1p = c1Rot.calculateEuclidianDistance(poiRot);

    double maxW = std::max(w1, std::max(w2, w3));
    double minW = std::min(w1, std::min(w2, w3));

    double wP = w1 + d1p / d1S * (w2 - w1 + d2S / d23 * (w3 - w2));

    if (wP < minW || wP > maxW) {
        DEBUG_STDOUT(
                "interpolated value " + std::to_string(wP) + " is out of bounds (" + std::to_string(minW) + "," +
                std::to_string(maxW) + ") "
                + std::to_string(d1p) + " " + std::to_string(d1S) + " " + std::to_string(d2S) + " " +
                std::to_string(d23));
    }

    return wP;
}

double Algorithms::interpolateBiquadraticOnTriangle(Coordinate3D pointOfInterest, Coordinate3D c1, double w1,
                                                    Coordinate3D gradC1, Coordinate3D c2, double w2,
                                                    Coordinate3D gradC2, Coordinate3D c3, double w3,
                                                    Coordinate3D gradC3) {
    //first step: rotate all points in the xy-layer
    Coordinate3D triVec1(c2 - c1);
    Coordinate3D triVec2(c3 - c1);

    Coordinate3D normVectorTriangle = triVec1.crossProduct(triVec2);
    Coordinate3D normVectorCoord = Coordinate3D{normVectorTriangle.x, normVectorTriangle.y, normVectorTriangle.z};
    SphericCoordinate3D normVectorSpheric = abm::util::toSphericCoordinates(normVectorCoord);
    double phiNorm = normVectorSpheric.phi;
    double thetaNorm = normVectorSpheric.theta;

    Coordinate3D poiRotPre = abm::util::rotateAboutThetaAndPhi(pointOfInterest, phiNorm, thetaNorm);
    Coordinate3D c1Rot = abm::util::rotateAboutThetaAndPhi(c1, phiNorm, thetaNorm);
    Coordinate3D poiRot{poiRotPre.x, poiRotPre.y, c1Rot.z};
    Coordinate3D c2Rot = abm::util::rotateAboutThetaAndPhi(c2, phiNorm, thetaNorm);
    Coordinate3D c3Rot = abm::util::rotateAboutThetaAndPhi(c3, phiNorm, thetaNorm);

    Coordinate3D gradC1Rot = abm::util::rotateAboutThetaAndPhi(gradC1, thetaNorm, phiNorm);
    Coordinate3D gradC2Rot = abm::util::rotateAboutThetaAndPhi(gradC2, thetaNorm, phiNorm);;
    Coordinate3D gradC3Rot = abm::util::rotateAboutThetaAndPhi(gradC3, thetaNorm, phiNorm);;

    //reorient vectors as a matter of curvature
    double gradC1Magn = gradC1Rot.getMagnitude();
    double gradC2Magn = gradC2Rot.getMagnitude();
    double gradC3Magn = gradC3Rot.getMagnitude();

    Coordinate3D gradC1RotRescale{gradC1Rot.x, gradC1Rot.y, 0};
    gradC1RotRescale.setMagnitude(gradC1Magn);
    Coordinate3D gradC2RotRescale{gradC2Rot.x, gradC2Rot.y, 0};
    gradC2RotRescale.setMagnitude(gradC2Magn);
    Coordinate3D gradC3RotRescale{gradC3Rot.x, gradC3Rot.y, 0};
    gradC3RotRescale.setMagnitude(gradC3Magn);


    //step two: find intersection of 2->3 and 1->poi
    double D = (c3Rot.y - c2Rot.y) * (poiRot.x - c1Rot.x) - (poiRot.y - c1Rot.y) * (c3Rot.x - c2Rot.x);

    double xS = ((c3Rot.x - c2Rot.x) * (poiRot.x * c1Rot.y - c1Rot.x * poiRot.y) -
                 (poiRot.x - c1Rot.x) * (c3Rot.x * c2Rot.y - c2Rot.x * c3Rot.y)) / D;
    double yS = ((c1Rot.y - poiRot.y) * (c3Rot.x * c2Rot.y - c2Rot.x * c3Rot.y) -
                 (c2Rot.y - c3Rot.y) * (poiRot.x * c1Rot.y - c1Rot.x * poiRot.y)) / D;
    Coordinate3D rS = Coordinate3D{xS, yS, c1Rot.z};

    //step three: interpolate the concentration quadratic
    double d23 = c2Rot.calculateEuclidianDistance(c3Rot);
    double d2S = c2Rot.calculateEuclidianDistance(rS);
    double d1S = c1Rot.calculateEuclidianDistance(rS);
    double d1p = c1Rot.calculateEuclidianDistance(poiRot);

    Coordinate3D d23_unit(c3Rot - c2Rot);
    d23_unit.setMagnitude(1.0);
    Coordinate3D d1S_unit(rS - c1Rot);
    d1S_unit.setMagnitude(1.0);

    double gradC1_1S = gradC1RotRescale.scalarProduct(d1S_unit);
    double gradC2_1S = gradC2RotRescale.scalarProduct(d1S_unit);
    double gradC3_1S = gradC3RotRescale.scalarProduct(d1S_unit);

    double gradCS_1S = gradC2_1S + d2S / d23 * (gradC3_1S - gradC2_1S);

    double a = (gradCS_1S - gradC1_1S) / (2.0 * d1S);
    double b = gradC1_1S;
    double c = w1;
    double wP = a * d1p * d1p + b * d1p + c;

    return wP;
}

using namespace boost::numeric::ublas;
bool Algorithms::InvertMatrix(const matrix<double> &input, matrix<double> &inverse) {
    typedef permutation_matrix<std::size_t> pmatrix;

    // create a working copy of the input
    matrix<double> A(input);

    // create a permutation matrix for the LU-factorization
    pmatrix pm(A.size1());

    // perform LU-factorization
    int res = lu_factorize(A, pm);
    if (res != 0)
        return false;

    // create identity matrix of "inverse"
    inverse.assign(identity_matrix<double>(A.size1()));

    // backsubstitute to get the inverse
    lu_substitute(A, pm, inverse);

    return true;
}

double Algorithms::meanAngle(double *angles, int size) {
    double y_part = 0, x_part = 0;
    int i;

    for (i = 0; i < size; i++) {
        x_part += cos(angles[i] * M_PI / 180);
        y_part += sin(angles[i] * M_PI / 180);
    }

    return atan2(y_part / size, x_part / size) * 180 / M_PI;
}