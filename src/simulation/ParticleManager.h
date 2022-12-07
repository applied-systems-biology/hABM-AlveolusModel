//  Copyright by Christoph Saffer, Paul Rudolph, Sandra Timme, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#ifndef PARTICLEMANAGER_H
#define PARTICLEMANAGER_H

#include <memory>
#include <vector>

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include "io/XMLFile.h"
#include "simulation/Particle.h"
#include "simulation/neighbourhood/StaticBalloonList.h"
#include "utils/io_util.h"

struct TRIANGLE3D {
    bool outside;
    unsigned int neighbourIds[3];
};

using namespace boost::numeric;

class ParticleManager {
public:
    /// Class for managing diffusion of particles inside the site
    explicit ParticleManager(Site *site);
    std::vector<std::shared_ptr<Particle>> &getAllParticles() { return allParticles; };

    /*!
     * Initializes particles (i.e. takes particles from fixed input file)
     * @param site Site object of environment (e.g. AlveoleSite)
     * @param parameters Parameter object that contains parameters of particle manager
     * @param input_dir String that contains input directory
     */
    void initializeParticles(Site *site,
                             const abm::util::SimulationParameters::ParticleManagerParameters &parameters,
                             const std::string &input_dir);

    /*!
     * Input of particles
     * @param time_delta Double for current timestep
     */
    void inputOfParticles(double time_delta);

    /*!
     * Creates a particle
     * @param site Site object of environment (e.g. AlveoleSite)
     * @param pos Coordinate3D object for position of particle
     * @param conc Double that contains initial concentration
     * @param dc Double that contains the diffusion coefficient
     * @return Particle object
     */
    std::shared_ptr<Particle> createParticle(Site *site, Coordinate3D pos, double conc, double dc);

    /*!
     * Replaces or removes a particle (i.e. for absorbing boundaries)
     * @param site
     * @param particle
     * @param newParticle
     */
    void replaceParticle(Site *site, Particle *particle, std::shared_ptr<Particle> newParticle);

    /// Inits clean up of all particles
    void cleanUpAllParticles();
    void computeMinTimestepDistribution();
    void diffusionPSE(double timestep, double time_delta);
    void includeParticleXMLTagToc(XMLFile *xmlTags);
    void setCleanChemotaxis(bool val) { clean_chemotaxis = val; };
    void setAECCells(std::vector<SphericCoordinate3D> AECT1, std::vector<SphericCoordinate3D> AECT2);

    std::shared_ptr<Particle> getParticleByPosition(Coordinate3D pos);
    StaticBalloonList *getParticleBalloonList() { return particleBalloonList.get(); }
    double getGradient(const Coordinate3D &position);
    double getSumChemokine();
    bool steadyStateReached(double current_time);
    [[nodiscard]] double getDiffusionCoefficient() const { return dc; };
    [[nodiscard]] bool getallowHigherDT() const { return allowHigherDT; };

private:
    void insertConcentrationAtArea(Site *site, double time_delta);
    void addTriangle(unsigned int id1, unsigned int id2, unsigned int id3);
    void triangulationFromDirectInput(Site *site,
                                      const abm::util::SimulationParameters::ParticleManagerParameters &parameters,
                                      const std::string &input_dir);
    void extractTriangles();
    int get_closest_AEC_ID(Coordinate3D position, int type);

    std::vector<std::shared_ptr<Particle>> allParticles;
    std::unique_ptr<StaticBalloonList> particleBalloonList;
    std::map<unsigned int, Particle *> particleById;
    unsigned int noOfParticles;
    unsigned int idHandling;
    ublas::compressed_matrix<double> P;
    ublas::vector<double> concBlas;
    ublas::vector<double> dConcBlas;
    std::shared_ptr<Particle> conidiumParticle;
    std::vector<std::shared_ptr<Particle>> aecParticles;
    std::vector<int> aecParticlesCells;
    std::vector<SphericCoordinate3D> alvEpithTypeOne;
    std::vector<SphericCoordinate3D> alvEpithTypeTwo;
    std::vector<TRIANGLE3D> triangles;
    std::string particleInputDelauneyFile;
    double **conc;
    double dc;
    double sumAreaAECParticles;
    double sumAreaAEcParticlesCells[100];
    double aecSecretionratePerGrid[100];
    double particleSecretionMoleculePerCellMin;
    bool allowHigherDT;
    bool clean_chemotaxis = true;
    bool drawIsolines;
    Site *site;

};

#endif    /* PARTICLEMANAGER_H */

