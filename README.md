# Hybrid Agent-Based Model (hABM) - Alveolus Model

An analytical and computational model framework for simulating early *Aspergillus fumigatus* infection scenarios in the human or murine alveolus.

Christoph Saffer, Paul Rudolph, Sandra Timme, Marco Blickensdorf, Johannes Pollmächer

Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge\
https://www.leibniz-hki.de/en/applied-systems-biology.html \
HKI-Center for Systems Biology of Infection\
Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Institute (HKI)\
Adolf-Reichwein-Straße 23, 07745 Jena, Germany

The project code is licensed under BSD 2-Clause.\
See the LICENSE file provided with the code for the full license.

## Project

Welcome to the repository of the Hybrid Agent-Based Model (hABM) - Alveolus Model developed by the **Applied Systems Biology** team of **Leibniz-HKI**.
The main purpose of this framework includes simulating an early *Aspergillus fumigatus* infection scenario in one human or murine alveolus. Moreover, we provide an analytical Surrogate Infection Model (SIM) that predicts hABM outcomes in a fraction of time and allows extended analyses of infection dynamics.

![](alveolusImage.png "To-scale representation of one murine (left) and human (right) Alveolus").

To-scale representation of the murine (left) and human (right) alveolus as a 3/4
sphere in the hABM. It is composed of an epithelial cell layer consisting of alveolar epithelial
cells (AEC) type I (yellow) and type II (blue). A single Aspergillus fumigatus conidium (red) is
randomly positioned on the inner alveolar surface and the contacting AEC is secreting
chemokines (white isolines) to attract alveolar macrophages (green) towards the conidium.
The alveolar entrance ring and the pores of Kohn (black) are the boundaries of the system.

## Requirements
- C++17 (or higher)
- Boost (>= Version 1.56)
- OpenMP (Version 4.5)
- PovRay (Version 3.7.0.8)
- CMake (recommended for build process)

## Getting started
### Build process

After cloning the repository to a local folder (here: abm/), the framework can be built via cmake:

`~/hABM-AlveolusModel$ mkdir build; cd build`

Release or Debug mode: 

`~/hABM-AlveolusModel/build$ cmake -DCMAKE_BUILD_TYPE=Release .. `
or `~/hABM-AlveolusModel/build$ cmake -DCMAKE_BUILD_TYPE=Debug .. `

`~/hABM-AlveolusModel/build$ make `

The compiled files can be found in the build/ folder.

### Run test configurations

To test if everything was compiled accordingly, run test configuration:

`~/hABM-AlveolusModel/build$ cd test/`

`~/hABM-AlveolusModel/build/test$ ./test_configurations`

(Test must be executed from the folder)

If the tests pass, the framework and its corresponding libraries were successfully installed.

### Model usage and input

Generally, the executable program `build/src/hABM` must be run with a `.json` configuration file as program argument

`~/hABM-AlveolusModel$ build$/src/hABM <config-file>.json` 

Additional input files 
- `analyser-config.json`
- `input-config.json`
- `output-config.json`
- `simulator-config.json`
- `visualisation-config.json`

are located in the folder specified as `config_path` variable in `<config-file>.json`.

The output is written to the folder specified in the `output_path` variable in the `<config-file>.json`. 

### Output formats

The output is written to `output_path/AlveolusOutput/SIMULATION_FOLDER/` and contains
- measurements (if activated in `analyser-config.json`)
- XML output (if activated in `output-config.json`)
- visualization (if activated in `visualization-config.json`)

The main output format for generating data applicable to analytical models is found in `measurements/` folder.

## Configurations

Predefined configurations including test scenarios can be found in the folder 'configurations/'. 

Human scenario: Config folder `configurations/namStudyHuman` and config file `config_human.json`
  
Murine scenario: Config folder `configurations/namStudyMouse` and config file `config_mouse.json` 

The parameters for each scenario can be adjusted in the `.json` files in the corresponding configuration folder/files.

## Run simulations 

To run simulations, execute:

Human: `~/hABM-AlveolusModel$ build/src/hABM config_human.json`

Mouse: `~/hABM-AlveolusModel$ build/src/hABM config_mouse.json`

For parameter screening, the cartesian product of all input sets in the `config_*.json` file in `parameter_screening` is calculated and the simulations are started for all parameter combinations.

## General structure
The framework is structured as followed:

- **Analysis_SurrogateInfectionModel**: Data processing and framework for analytical Surrogate Infection Model (SIM)
- **configurations**: configuration folders including `.json`files that define single scenarios
- **input**: precalculated data or parameters for the model
- **src**: source `.cpp` and header `.h` files of the framework
- **test**: test scenarios and unit tests using doctest library
- **cmake**: CMake specific files

For a detailed overview, we refer to the doxygen documentation (Doxygen-Output.zip).

## Previous publications

The essence of the hybrid agent-based model was applied to various studies about *Aspergillus fumigatus* lung infection in the past, such as:
- [Pollmächer *et. al.* (2014)](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0111630)
- [Pollmächer *et. al.* (2015)](https://www.frontiersin.org/articles/10.3389/fmicb.2015.00503/)
- [Blickensdorf *et. al.* (2019)](https://www.frontiersin.org/articles/10.3389/fimmu.2019.00142/)
- [Blickensdorf *et. al.* (2020)](https://www.frontiersin.org/articles/10.3389/fmicb.2020.01951/)
