//  Copyright by Christoph Saffer, Paul Rudolph, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#include "simulation/NeighbourhoodLocator.h"
#include "simulation/Agent.h"
#include "simulation/neighbourhood/Collision.h"


NeighbourhoodLocator::NeighbourhoodLocator(Site *site) {
    site_ = site;
}

NeighbourhoodLocator::~NeighbourhoodLocator() {
}

void NeighbourhoodLocator::instantiate() {
}

std::vector<std::shared_ptr<Collision>> NeighbourhoodLocator::getCollisions(Agent *agent) {
    return {};
}

void NeighbourhoodLocator::updateDataStructures(SphereRepresentation *sphereRep) {
}

void NeighbourhoodLocator::setThresholdDistance(double threshDistance) {
    thresholdDistance = threshDistance;
}

void NeighbourhoodLocator::removeSphereRepresentation(SphereRepresentation *sphereRep) {
}

void NeighbourhoodLocator::addSphereRepresentation(SphereRepresentation *sphereRep) {
}

std::string NeighbourhoodLocator::getTypeName() {
    return "NeighbourhoodLocator";
}

std::vector<Coordinate3D>
NeighbourhoodLocator::getCollisionSpheres(SphereRepresentation *sphereRep, Coordinate3D dirVec) {
    return std::vector<Coordinate3D>();
}

