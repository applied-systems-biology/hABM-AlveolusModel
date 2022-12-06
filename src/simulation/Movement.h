//  Copyright by Christoph Saffer, Paul Rudolph, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#ifndef ABM_MOVEMENT_H_
#define    ABM_MOVEMENT_H_

#include <memory>

#include "basic/Coordinate3D.h"
#include "basic/Randomizer.h"


class Site;
class Movement {
public:
  // Abstract class for movement action of cells
    explicit Movement(unsigned int spatial_dimensions);
    virtual ~Movement() = default;

    void setSite(Site *site);
    void setCurrentPosition(Coordinate3D *agent_pos);
    void setCurrentMove(Coordinate3D *);
    void setCurrentTimestep(double timestep);
    virtual void setPreviousMove(Coordinate3D *) {}
    [[nodiscard]] Coordinate3D *getCurrentMove();
    [[nodiscard]] double getCurrentTimestep() const { return current_timestep_; }
    virtual double getSpeed() { return 0; }
    virtual double getStartingTime();
    virtual std::string getMovementName();
    virtual Coordinate3D *move(double, double diffusion_constant);
    virtual void annulatePersistence() {}

protected:
    unsigned int spatial_dimensions_;
    Site *site_{};
    Coordinate3D *current_pos_{};
    std::unique_ptr<Coordinate3D> current_move_{};
    double current_timestep_{};

};

#endif    /* ABM_MOVEMENT_H_ */

