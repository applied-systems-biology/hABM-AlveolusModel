//  Copyright by Christoph Saffer, Paul Rudolph, Sandra Timme, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#include "Ingestion.h"
#include "simulation/Interaction.h"
#include "simulation/Cell.h"
#include "simulation/Site.h"
#include "simulation/Interactions.h"


void Ingestion::handleInteraction(Interaction *interaction, Cell *cell, double timestep, double current_time) {
    //active cell is the ingested one (the one that is eaten)
    //passive cell is the eater-cell

    Cell *passiveCell, *activeCell;

    if (isIngestingType(cell)) {
        activeCell = interaction->getOtherCell(cell);
        passiveCell = cell;
    } else {
        activeCell = cell;
        passiveCell = interaction->getOtherCell(cell);
    }

    Coordinate3D shift = passiveCell->getEffectiveConnection(activeCell);
    double currentDistance = shift.getMagnitude();

    // put the ingested cell at r/2 distance of the center of the ingesting cell
    double multiplicatorForNormalVec =
            passiveCell->getSurface()->getAllSpheresOfThis().front()->getRadius() / (2 * currentDistance) - 1;
    shift *= multiplicatorForNormalVec;
    activeCell->shiftPosition(&shift, current_time, 0, getTypeName());

    // prevent that ingested cell goes out of the site
    Site *currentSite = cell->getSite();
    std::string siteName = currentSite->getType();

    // Ensure that passive cell does not get out of spheric site
    double r = currentSite->getRadius();
    double border = currentSite->getThicknessOfBorder();
    SphericCoordinate3D currentPos = abm::util::toSphericCoordinates(activeCell->getPosition());
    if (currentPos.r < r - border / 2) {
        currentPos.r = 0.5 * (abm::util::toSphericCoordinates(passiveCell->getPosition()).r + r - border / 2);
        activeCell->setPosition(abm::util::toCartesianCoordinates(currentPos));
    } else if (currentPos.r > r + border / 2) {
        currentPos.r = 0.5 * (abm::util::toSphericCoordinates(passiveCell->getPosition()).r + r + border / 2);
        activeCell->setPosition(abm::util::toCartesianCoordinates(currentPos));
    }
    if (activeCell->getCurrentCellState()->getStateName() != "Death") {
        // Conidia is taken up alive, then it is removed from the list which induces the particle cleanup
        currentSite->getAgentManager()->removeConidiaFromList(activeCell->getId(), current_time);
    }


    std::shared_ptr<Collision> currentCollision;
    while ((currentCollision = interaction->getNextCollision()) != 0) {
        currentCollision = nullptr;
    }
}

bool Ingestion::isIngestingType(Cell *cell) {
    bool result = false;

    if (cell->getTypeName().compare("Macrophage") == 0) { result = true; }

    return result;
}

std::string Ingestion::getTypeName() const {
    return "Ingestion";
}