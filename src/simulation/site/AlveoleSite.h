//  Copyright by Christoph Saffer, Paul Rudolph, Sandra Timme, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#ifndef ALVEOLESITE_H
#define    ALVEOLESITE_H

#include <vector>
#include <utility>
#include <algorithm>

#include "SphereSite.h"
#include "basic/SphericCoordinate3D.h"


class AlveoleSite : public SphereSite {
public:
    /// Class that constructs the alveolus environment in the human and murine case
    AlveoleSite(abm::util::SimulationParameters::SiteParameters *parameters,
                double time_delta,
                unsigned int spatial_dimensions,
                const std::string &input_dir,
                Randomizer *random_generator, std::shared_ptr<InSituMeasurements> measurements);

    /*!
     * Generates random direction vector on alveoleSite (sphere)
     * @param position Coordinate3D object that contains current position
     * @param length Double that containts desired length of random direction
     * @return Coordinate3D object that contains direction vector
     */
    Coordinate3D generateRandomDirectionVector(Coordinate3D position, double length) final;

    /// Same as generateRandomDirectionVector but based on previous movement
    Coordinate3D generatePersistentDirectionVector(Coordinate3D position, double length, Coordinate3D prevVector, double previousAlpha) final;

    /// Same as generateRandomDirectionVector but with a random bias (e.g. molecules)
    Coordinate3D generateBiasedRandomDirectionVector(Agent *agent, Coordinate3D position, double length) final;

    bool overPOK(Coordinate3D position);
    bool onAECTObstacleCell(Coordinate3D position) final;
    bool containsPosition(Coordinate3D position) final;
    bool overAECT1(SphericCoordinate3D posConida) final;
    void includeSiteXMLTagToc(XMLFile *xmlFile) const final;
    double getFeatureValueByName(std::string name) final;
    double getThicknessOfBorder() final { return thicknessOfBorder; }
    double getLowerThetaBound() final { return thetaLowerBound; }
    double getDistanceFromBoundary(Coordinate3D position) final;
    Coordinate3D getRandomPosition() final;
    Coordinate3D getRandomBoundaryPoint() final;
    Coordinate3D getRandomMinDistanceToBoundaryPosition(double minDistanceToBoundary) final;
    [[nodiscard]] std::string getType() const final { return "AlveoleSite"; }
    std::vector<SphericCoordinate3D> getAECT1() final { return alvEpithTypeOne; };
    std::vector<SphericCoordinate3D> getAECT2() final { return alvEpithTypeTwo; };
    std::vector<SphericCoordinate3D> getPOK() final { return poresOfKohn; };
    friend void InSituMeasurements::observeMeasurements(const SimulationTime &time);

private:
    void initializeAgents(const abm::util::SimulationParameters::AgentManagerParameters &parameters,
                          const std::string &input_dir,
                          double current_time, double time_delta) override;

    /// Inserts the alveolar epithelium type 1 cells in the system
    void includeAlveolarEpitheliumType1();
    /// Inserts the pores of kohn in the system in the human case
    void includePoresOfKohnHuman();
    /// Inserts the  alveolar epithelium type 2 in the system in the human case
    void includeAlveolarEpitheliumType2Human();
    /// Inserts the pores of kohn in the system in the murine case
    void includePoresOfKohnMouse();
    /// Inserts the  alveolar epithelium type 2 in the system in the murine case
    void includeAlveolarEpitheliumType2Mouse();

    Coordinate3D generateDirectedVector(Coordinate3D position, SphericCoordinate3D posOfGoal, double length) final;
    Coordinate3D generateDirectedVector(Coordinate3D position, double alpha, double length) final;
    static double retrieveDirectionAngleAlpha(SphericCoordinate3D ownPos, SphericCoordinate3D goalPos);
    double minDistanceToPoK(const SphericCoordinate3D &sc3d);
    void calculateCrossPoints();

    int organism{};
    bool respirationEnabled{};
    bool obstacleIsOnType1{};
    int noOfPoK{};
    int noOfAEC2{};
    double thicknessOfBorder{};
    double surfactantThickness{};
    double thetaLowerBound{};
    double siteRadiusMin{};
    double radiusAlvEpithTypeOne{};
    double lengthAlvEpithTypeTwo{};
    double radiusPoresOfKohn{};
    double r0AEC1{};
    double r0AEC2{};
    double opR{};

    SphericCoordinate3D posObstacle{};
    SphericCoordinate3D cellOfObstacle{};
    std::vector<std::pair<int, int> > pairCells{};
    std::vector<SphericCoordinate3D> alvEpithTypeOne{};
    std::vector<SphericCoordinate3D> alvEpithTypeTwo{};
    std::vector<SphericCoordinate3D> poresOfKohn{};
    std::vector<Coordinate3D> crossPoints{};
};

#endif    /* ALVEOLESITE_H */

