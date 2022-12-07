//  Copyright by Christoph Saffer, Paul Rudolph, Sandra Timme, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#include "simulation/Condition.h"
#include "simulation/Cell.h"


Condition::Condition() {
}

Condition::Condition(Cell *cell) {
    this->cell = cell;
    condition = "";
}

Condition::Condition(std::string condition) {
    this->condition = condition;
    cell = 0;
}

Condition::Condition(const Condition &orig) {
}

Condition::~Condition() {
}

Cell *Condition::getCell() {
    return cell;
}

std::string Condition::getStringCondition() {
    return condition;
}

bool Condition::isFulfilled(Condition *condition) {
    bool fulfilled = false;
    if (cell != 0) {
        if (condition->getCell() == 0) {
            fulfilled = (condition->getStringCondition().compare(cell->getTypeName()) == 0);
        } else {
            fulfilled = (cell == condition->getCell());
        }

    } else {
        if (condition->getCell() == 0) {
            fulfilled = (condition->getStringCondition().compare(this->condition) == 0);
        } else {
            fulfilled = (condition->getCell()->getTypeName().compare(this->condition) == 0);
        }
    }
    return fulfilled;
}