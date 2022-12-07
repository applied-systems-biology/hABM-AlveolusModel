//  Copyright by Christoph Saffer, Paul Rudolph, Sandra Timme, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#ifndef SAMPLER_H
#define    SAMPLER_H

#include <iostream>
#include <vector>

#include "basic/Randomizer.h"

class Sampler {
public:
  // Class for randomly sampling values
    Sampler();

    Sampler(const char *);
    virtual ~Sampler() = default;
    virtual void setSampleFunction(int);
    virtual void setSampleRange(double, double, unsigned int dim = 1);
    virtual void sample(Randomizer *randomizer, unsigned int distrColumn = 0);
    virtual double getSampledValue(unsigned int dim = 1);

protected:
    int samplingBasis;
    const char *filename;
    int samplingFunction;
    std::vector<double> lowerBound;
    std::vector<double> upperBound;

};

#endif    /* SAMPLER_H */

