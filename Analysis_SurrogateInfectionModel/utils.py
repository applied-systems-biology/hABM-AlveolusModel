#   Copyright by Christoph Saffer
#   Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
#   https://www.leibniz-hki.de/en/applied-systems-biology.html
#   HKI-Center for Systems Biology of Infection
#   Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
#   Adolf-Reichwein-Straße 23, 07745 Jena, Germany
#
#   This code is licensed under BSD 2-Clause
#   See the LICENSE file provided with this code for the full license.

import pandas as pd
import numpy as np
import random
from scipy.optimize import minimize
import warnings
warnings.filterwarnings('ignore')  # setting ignore as a parameter


# Process infection scores and input variables to provide correct input format
def retrieve_infection_scores_and_input_variables(file_name: str, cond_Dc=None, cond_sAEC=None, cond_nOfCon=None) -> list:
    df = pd.read_csv(file_name, index_col=0)

    # Input variables
    nOfMs = np.sort(df["nOfM"].unique())

    if cond_Dc is None:
        Dcs = np.sort(df["Dc"].unique())
        rwdataDC = [6000]  # Quasi-random walk
    else:
        Dcs = [cond_Dc]
        rwdataDC = []

    if cond_sAEC is None:
        sAECs = np.sort(df["sAEC"].unique())[2:8]
        rwdatasAEC = [1, 10]  # Quasi-random walk
    else:
        sAECs = [cond_sAEC]
        rwdatasAEC = []

    if cond_nOfCon is None:
        nOfCons = np.sort(df["nOfCon"].unique())
    else:
        nOfCons = [cond_nOfCon]

    infScores = []
    variables = []

    for nOfCon in nOfCons:
        for sAEC in sAECs:
            for Dc in Dcs:
                for nOfM in nOfMs:
                    infScore = df[df["nOfCon"] == nOfCon][df["Dc"] == Dc][df["sAEC"] == sAEC][df["nOfM"] == nOfM]["infScore"].values[0]
                    infScores.append(infScore)
                    variables.append([nOfM, nOfCon, sAEC, Dc])

        for Dc in rwdataDC:
            for sAEC in rwdatasAEC:
                for nOfM in nOfMs:
                    infScore = df[df["nOfCon"] == nOfCon][df["Dc"] == Dc][df["sAEC"] == sAEC][df["nOfM"] == nOfM]["infScore"].values[0]
                    infScores.append(infScore)
                    variables.append([nOfM, nOfCon, sAEC, Dc])

    infScores = np.array(infScores)
    variables = np.array(variables).transpose()

    return [infScores, variables]


# Retrieve seed points to optimize SIM derived from generalized functions
def retrieve_seed_point_for_SIM(generalized_fun, system: str, calculate_seeds_from_cef_para=False) -> list:

    # read in cef parameters
    cef_parameter = pd.read_csv("fitted_parameter/CEF_"+system+".csv", index_col=0)
    # seed points for cef optimization
    seedpoint_beta = np.array([1, -15, -5, 1, 1, -1, -1])  # according to initial analysis of generalized functions
    seedpoint_gamma = np.array([1, -5, -15, 1, 1, 1, 1])  # according to initial analysis of generalized functions

    # prepare data for optimization
    betas = cef_parameter["beta"]
    gammas = cef_parameter["gamma"]
    input_variables = np.array(cef_parameter[['nOfCon', 'sAEC', 'Dc']]).transpose()

    def objective_function_betas(parameters: list) -> float:
        list_of_squared_errors = ((betas - generalized_fun(input_variables, parameters)) ** 2)
        return np.sum(list_of_squared_errors) / len(list_of_squared_errors)

    def objective_function_gammas(parameters: list) -> float:
        list_of_squared_errors = ((gammas - generalized_fun(input_variables, parameters)) ** 2)
        return np.sum(list_of_squared_errors) / len(list_of_squared_errors)

    def optimize_generalized_functions(objective_function, seed, number_of_trials, noise_to_seed):
        optimal_value = objective_function(seed)
        optimal_parameters = seed
        for i in range(number_of_trials):
            current_seed = seed + np.array([random.gauss(0, noise_to_seed) for p in range(0, len(seed))])
            result = minimize(objective_function, current_seed, method='tnc')
            if result.fun < optimal_value:
                optimal_value = result.fun
                optimal_parameters = result.x

        return optimal_parameters

    if calculate_seeds_from_cef_para:
        # optimize beta parameters
        beta_opt = optimize_generalized_functions(objective_function_betas, seedpoint_beta, number_of_trials=50, noise_to_seed=0.1)
        # optimize gamma parameters
        gamma_opt = optimize_generalized_functions(objective_function_gammas, seedpoint_gamma, number_of_trials=50, noise_to_seed=0.1)
        seed_for_sim = []
        # combine beta and gamma parameters to SIM seed parameters
        for x in beta_opt:
            seed_for_sim.append(x)
        for x in gamma_opt:
            seed_for_sim.append(x)

        print("Found seedpoint for SIM by optimizing generalized functions.")

    else:
        # Mainly used for cross validation
        # Optimal seed points derived from our previous analysis of generalized functions fitted to CEF parameters:
        if system == "human":
            seed_for_sim = np.array([0.98067888, -14.71456769, -5.30023224, 1.29843273,
                                     0.06061795, -1.25441528, -3.35566329,
                                     0.15103859, -2.04186545, -10.54154265, 1.01549763,
                                     1.05436988, 1.29116851, 0.54643654])

        else:
            seed_for_sim = np.array([7.4686941, -12.69401084, -4.33377718, 1.11667974,
                                     1.32255414, -0.15714425, -0.68052268,
                                     0.18548205, -4.44837653, -19.60556799, 1.79473858,
                                     1.07459942, 0.92629279, 0.35437558])

    return seed_for_sim


# Read in optimal SIM parameters
def read_parameters(path: str) -> list:
    df = pd.read_csv(path, index_col=0)
    p = np.array(df.iloc[0, :])
    # Due to numerical reasons, we perform a variable transformation in the logistic functions
    # see generalized_function_f in train_surrogate_infection_model.py
    # (x2or3 s/D)^x4 = exp(x4 * log(s/D) + x4 * log(x2or3)) ==> x2or3' = x4 * log(x2or3)
    p[1] = p[3] * np.log(p[1])
    p[2] = p[3] * np.log(p[2])
    p[8] = p[10] * np.log(p[8])
    p[9] = p[10] * np.log(p[9])
    return p


# Save optimal SIM parameters
def write_parameters(path: str, p: np.array):
    # Due to numerical reasons, we perform a variable transformation in the logistic functions
    # see generalized_function_f in train_surrogate_infection_model.py
    # exp(x4 * np.log(s/D) + x2or3') = (exp(x2or3')^(1/x4) * s/D)^(x4) ==> x2or3 = exp(x2or3')^(1/x4)
    p[1] = np.exp(p[1])**(1/p[3])
    p[2] = np.exp(p[2])**(1/p[3])
    p[8] = np.exp(p[8])**(1/p[10])
    p[9] = np.exp(p[9])**(1/p[10])

    df = pd.DataFrame([p], columns=["b1", "b2", "b3", "b4", "b5", "b6", "b7", "g1", "g2", "g3", "g4", "g5", "g6", "g7"])
    df.to_csv(path)
