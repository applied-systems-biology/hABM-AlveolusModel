//  Copyright by Christoph Saffer, Paul Rudolph, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#ifndef CONDITION_H
#define	CONDITION_H

#include <string>

class Cell;

class Condition {
public:
  // Class for conditioning interactions or states to environmental or cell-specific factors.
    Condition();
    
    explicit Condition(Cell* cell);
    explicit Condition(std::string condition);
    Condition(const Condition& orig);
    virtual ~Condition();
    bool isFulfilled(Condition* condition);
    Cell* getCell();
    std::string getStringCondition();

private:
    Cell* cell;
    std::string condition;
    
};

#endif	/* CONDITION_H */

