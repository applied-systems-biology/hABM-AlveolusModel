#   Copyright by Christoph Saffer
#   Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
#   https://www.leibniz-hki.de/en/applied-systems-biology.html
#   HKI-Center for Systems Biology of Infection
#   Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
#   Adolf-Reichwein-Straße 23, 07745 Jena, Germany
#
#   This code is licensed under BSD 2-Clause
#   See the LICENSE file provided with this code for the full license.

from scipy.optimize import minimize

from utils import *


# Global variables
optimize = False  # if False, use stored values to calculate MAD # if True, optimizer is started
save_to_file = False  # save parameters to file


# Generalized functions f for f_beta and f_gamma
def generalized_function_f(variables: list, parameters: list) -> float:
    nOfCon, sAEC, Dc = variables
    ratio = np.log(sAEC / Dc)
    x1, x2, x3, x4, x5, x6, x7 = parameters

    # Ensure x1 and x4 to be positive
    x1 = np.sqrt(x1 ** 2)
    x4 = np.sqrt(x4 ** 2)

    # Due to numerical reasons, we perform a variable transformation in the logistic functions
    # For x2 = exp(x2')^(1/x4), x3 = exp(x3')^(1/x4), we retrieve the original formula
    logistic_functions = (1 / (1 + np.exp(x4 * ratio + x2))) - (1 / (1 + np.exp(x4 * ratio + x3)))
    # To provide the correct parameter values, the variables are transformed in the read/write functions in utils.py

    return x1 * logistic_functions * nOfCon.astype(float) ** x6 + x5 * nOfCon.astype(float) ** x7


# Surrogate infection model (SIM)
def surrogate_infection_model(variables: list, parameters: list) -> float:
    nOfM, nOfCon, sAEC, Dc = variables
    beta = parameters[:7]
    gamma = parameters[7:]
    gf_input = [nOfCon, sAEC, Dc]

    exponent = - generalized_function_f(gf_input, beta) * nOfM ** generalized_function_f(gf_input, gamma)

    return np.exp(exponent)


# Objective function to calibrate SIM parameters
def objective_function_to_minimize(parameters: np.array) -> float:
    list_of_squared_errors = ((infection_scores - surrogate_infection_model(input_variables, parameters)) ** 2)
    return np.sum(list_of_squared_errors)/len(list_of_squared_errors)


# Calculate mean average deviation (MAD) for SIM
def calculate_mad(parameters: np.array) -> float:
    list_of_absolute_errors = np.abs(infection_scores - surrogate_infection_model(input_variables, parameters))
    return np.sum(list_of_absolute_errors)/len(list_of_absolute_errors)


# Optimizer for SIM using Scipy.optimize.minimze package
def train_surrogate_infection_model(seed: np.array, number_of_trials: int, noise_to_seed: float, number_of_second_trails=10) -> np.array:
    optimal_value = objective_function_to_minimize(seed)
    optimal_parameters = seed
    for j in range(number_of_second_trails):
        # Change seedpoint for next runs to previous optimum
        seed = optimal_parameters
        for i in range(number_of_trials):
            # Add gaussian noise to the current seedpoint to potentially find better optimum
            current_seed = seed + np.array([random.gauss(0, noise_to_seed) for p in range(0, len(seed))])
            # Optimize function
            result = minimize(objective_function_to_minimize, current_seed, method='tnc')
            if result.fun < optimal_value:
                optimal_value = result.fun
                optimal_parameters = result.x
                print("Run:", j*number_of_trials+(i+1), "/", number_of_second_trails*number_of_trials, ": Found new optimum for LSE =", optimal_value)

    return optimal_parameters, optimal_value


if __name__ == '__main__':

    # Optimize for both systems
    for system in ["human", "mouse"]:
        # Read in processed data
        input_file = "processed_data/" + system + "_processed_data.csv"  # path to retrieved infection scores

        # Retrieve infection scores and input variables from processed data
        infection_scores, input_variables = retrieve_infection_scores_and_input_variables(input_file)
        print("\n#### Optimize", system, "system:\nNumber of datapoints: ", len(infection_scores))

        if optimize:
            # Retrieve seedpoints derived from generalized functions fitted to CEF parameters
            parameters = np.array(retrieve_seed_point_for_SIM(generalized_function_f, system, calculate_seeds_from_cef_para=True))
            print("LSE from seedpoint: ", objective_function_to_minimize(parameters))
            print("MAD from seedpoint: ", calculate_mad(parameters))
            # Train surrogate model
            p, x = train_surrogate_infection_model(parameters, number_of_trials=100, noise_to_seed=0.5)
        else:
            p = read_parameters("fitted_parameter/SIM_"+system+".csv")

        print("SIM", system, "MAD =", calculate_mad(p))

        # Save parameters to file
        if save_to_file:
            write_parameters("fitted_parameter/SIM_"+system+".csv", p)

        # EXAMPLE: Usage of SIM for datapoint [nOfM=2.0, nOfCon=1, sAEC=1500, Dc=20] and optimal parameters 'parameters':
        # infectionScore = surrogate_infection_model([2.0, 1, 1500, 20], optimal_parameters)