//  Copyright by Christoph Saffer, Paul Rudolph, Sandra Timme, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#include "basic/InversionSampler.h"

InversionSampler::InversionSampler() : Sampler() {
}

InversionSampler::InversionSampler(const char *filename) : Sampler(filename) {

}

InversionSampler::InversionSampler(const InversionSampler &orig) : Sampler(orig) {
}

InversionSampler::~InversionSampler() {
}

double InversionSampler::sampleFunction(double random) {
    double output;
    switch (samplingFunction) {
        case 0: //0.5*sin(x) in range [0,PI]
            setSampleRange(0.0, 3.141528);
            output = acos(1 - 2 * random);
            break;
        default: //uniformly distributed
            output = random;
            break;
    }
    return output;
}

void InversionSampler::sample(Randomizer *randomizer, unsigned int distrColumn) {
    sampledValue = sampleFunction(randomizer->generateDouble(1.0));
}

double InversionSampler::getSampledValue(unsigned int dim) {
    return sampledValue;
}