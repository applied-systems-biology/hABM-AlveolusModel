//  Copyright by Christoph Saffer, Paul Rudolph, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#ifndef BIASEDPERSISTENTRANDOMWALK_H
#define    BIASEDPERSISTENTRANDOMWALK_H

#include <math.h>

#include "simulation/Movement.h"
#include "basic/Randomizer.h"
#include "basic/Sampler.h"


class Agent;
class BiasedPersistentRandomWalk : public Movement {
public:
  // Class for the biased random walk with a persistent time. The bias is given by the cumulative gradient variable in the macrophage class.
    BiasedPersistentRandomWalk(Agent *agent, double persistenceTime, double speed, unsigned int spatial_dimensions);
    Coordinate3D *move(double, double dc) final;
    double getStartingTime() final;
    void setPreviousMove(Coordinate3D *) final;
    void annulatePersistence() final { persistence_time_left_ = 0; };
    void decrementLeftTime(double);
    void setNewPersistence();
    [[nodiscard]] bool persistentMove() const;
    Coordinate3D *movePersistent(double);
    Coordinate3D *moveBiasedRandomly(double);

private:
    Agent *agent_{};
    double persistence_time_{};
    double persistence_time_left_{};
    double persistence_time_start_{};
    double speed_{};
    double persistent_angle_alpha_2_d_{};
    std::unique_ptr<Coordinate3D> persistence_direction_{};
    std::unique_ptr<Sampler> sampler_{};
    std::shared_ptr<Coordinate3D> current_velocity_{};
};

#endif    /* BIASEDPERSISTENTRANDOMWALK_H */

