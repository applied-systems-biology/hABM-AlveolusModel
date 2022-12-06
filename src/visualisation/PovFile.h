//  Copyright by Christoph Saffer, Paul Rudolph, Marco Blickensdorf, Johannes Pollmächer
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#ifndef POVFILE_H
#define    POVFILE_H

#include <string>
#include <vector>

class Site;

class PovFile {
public:
    explicit PovFile(const std::string &filebody);

    void setGlobalPart(std::string);
    void setCameraPart(std::string);
    void setLightPart(std::string);
    void addLightPart(const std::string &);
    void setDimensions(std::string px_width, std::string px_height);
    void setBorderPart(std::string);
    void setBackgroundPart(std::string);
    void addPovObject(const std::string &);
    void doPovProcess(double current_time, std::string zpos = "");
    void transcribeSite(const Site &site);
    void transcribeParticles(const Site &site);
    void transcribeAgents(const Site &site);
    void transcribeAgents(const Site &site, double z_start, double z_end);

private:
    std::string global_part_;
    std::string camera_part_;
    std::string light_part_;
    std::string px_width_;
    std::string px_height_;
    std::string border_part_;
    std::string background_part_;
    std::string filebody_pov_;
    std::vector<std::string> pov_objects_;
};

#endif    /* POVFILE_H */

