//  Copyright by Christoph Saffer, Paul Rudolph, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#include "basic/Sampler.h"

Sampler::Sampler() {
    samplingBasis = 1;
    lowerBound.resize(3, 0);
    upperBound.resize(3, 0);

}

Sampler::Sampler(const char *filename) {
    this->filename = filename;
    samplingBasis = 0;
    //do normalisation and interpolation of the input data etc.
    lowerBound.resize(3, 0);
    upperBound.resize(3, 0);
}

void Sampler::setSampleFunction(int fct) {
    samplingFunction = fct;
    //std::cout << "sample fct set" << std::endl;
}

void Sampler::setSampleRange(double lowerBound, double upperBound, unsigned int dim) {
    this->lowerBound.at(dim - 1) = lowerBound;
    this->upperBound.at(dim - 1) = upperBound;
}

void Sampler::sample(Randomizer *randomizer, unsigned int distrColumn) {

}

double Sampler::getSampledValue(unsigned int dim) {
    return 0;
}