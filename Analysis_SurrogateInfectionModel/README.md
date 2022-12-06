## Data processing of simulation data

### Executable scripts

To reproduce results, run the following scripts with python 3.8. For each script, global variables at the top of the files have to be adjusted to produce desired outputs.

- `process_simulation_output.py` - processes the output data from the hABM and makes it accessible for analytical models


- `train_weibull_surival_model.py` - trains the WSM

- `train_compressed_exponential_function.py` - trains the CEF

- `train_surrogate_infection_model.py` - trains the SIM


- `run_cross_validation.py` - runs the 5 times 6-fold cross validation comparing SIM, MLP and RDF

- `run_symmetry_error_calculation.py` - calculates the symmetry error for identical ratios sAEC/Dc

### Additional classes and functions

- `utils.py` - contains functions that are used during the training of the models

- `MultilayerPerceptron.py` - contains the class for the Multilayer Perceptron (MLP)

### Folders

- `processed_data/` - contains the output from `process_simulation_output.py` and includes derived Infection Scores (IS)

- `fitted_parameters/` - contains the parameters from the models, e.g. from the `train_*.py` scripts 

- `cross_validation/` - contains the output from `run_cross_validation.py`

### Needed python packages (version)

- numpy (1.22.0)
- pandas (1.3.5)
- scikit-learn (1.0.2) 
- torch (1.10.1)
