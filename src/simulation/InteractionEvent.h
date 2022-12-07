//  Copyright by Christoph Saffer, Paul Rudolph, Sandra Timme, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#ifndef INTERACTIONEVENT_H
#define    INTERACTIONEVENT_H

#include <string>

class Interaction;

class InteractionEvent {
public:
  // Class for handling for default actions an interaction undertakes after being triggered.
    InteractionEvent(std::string previousState, std::string nextState);
    InteractionEvent(std::string previousState, std::string nextState, Interaction *interaction);

    void setDescriptiveName();
    std::string getDescriptiveName();
    std::string getNextState();
    std::string getPreviousState();
    Interaction *getInteraction();

private:
    std::string previousState;
    std::string nextState;
    std::string descriptiveName;
    Interaction *interaction;

};

#endif    /* INTERACTIONEVENT_H */

