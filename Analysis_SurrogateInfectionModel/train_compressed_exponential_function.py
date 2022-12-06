#   Copyright by Christoph Saffer
#   Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
#   https://www.leibniz-hki.de/en/applied-systems-biology.html
#   HKI-Center for Systems Biology of Infection
#   Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
#   Adolf-Reichwein-Straße 23, 07745 Jena, Germany
#
#   This code is licensed under BSD 2-Clause
#   See the LICENSE file provided with this code for the full license.

import numpy as np
import pandas as pd
from scipy.optimize import minimize

from utils import retrieve_infection_scores_and_input_variables


# Global variables
optimize = False  # if False, use stored values to calculate MAD # if True, optimizer is started
save_to_file = False  # save parameters to file


# Compressed exponential function (CEF)
def compressed_exponential_function(variables: list, parameters: list) -> float:
    nOfM, nOfCon, sAEC, Dc = variables
    beta, gamma = parameters
    exponent = - beta * nOfM ** gamma

    return np.exp(exponent)


# Objective function to calibrate CEF parameters
def objective_function_to_minimize(parameters: np.array) -> float:
    list_of_squared_errors = ((infection_scores - compressed_exponential_function(input_variables, parameters)) ** 2)
    return np.sum(list_of_squared_errors)/len(list_of_squared_errors)


# Calculate mean average deviation (MAD) for CEF
def calculate_mad(parameters: np.array) -> float:
    list_of_absolute_errors = np.abs(infection_scores - compressed_exponential_function(input_variables, parameters))
    return np.sum(list_of_absolute_errors)/len(list_of_absolute_errors)


if __name__ == '__main__':

    for system in ["human", "mouse"]:
        # Set parameters for the corresponding system
        if system == "human":
            nOfMs = np.arange(2, 52, 2)  # human number of AM
            nOfCons = [1, 2]  # including high fungal burden: human number of conidia
        else:
            nOfMs = np.arange(0.1, 2.6, 0.1)  # mouse number of AM
            nOfCons = [1, 2, 3]  # including high fungal burden: mouse number of conidia
        sAECs = [1500, 5000, 15000, 50000, 150000, 500000]  # secretion rates
        Dcs = [20, 60, 200, 600, 2000, 6000]  # diffusion coefficients

        # Read in processed data
        input_file = "processed_data/" + system + "_processed_data.csv"  # path to retrieved infection scores
        res, counter, seed = 0, 0, [1, 1]
        fitted_parameter = []

        # Read optimal parameters
        df = pd.read_csv("fitted_parameter/CEF_"+system+".csv", index_col=0)

        # Loop over all nOfCon, Dc, and sAEC and calibrate CEF
        for nOfCon in nOfCons:
            for Dc in Dcs:
                for sAEC in sAECs:
                    infection_scores, input_variables = retrieve_infection_scores_and_input_variables(input_file, cond_Dc=Dc, cond_sAEC=sAEC, cond_nOfCon=nOfCon)

                    if optimize:
                        result = minimize(objective_function_to_minimize, seed).x
                        print("nOfCon =", nOfCon, ", Dc =", Dc, ", sAEC =", sAEC, ":: Fitted parameter beta, gamma =", result[0], ",", result[1])
                    else:
                        result = [df[df["nOfCon"] == nOfCon][df["Dc"] == Dc][df["sAEC"] == sAEC]["beta"].values[0], df[df["nOfCon"] == nOfCon][df["Dc"] == Dc][df["sAEC"] == sAEC]["gamma"].values[0]]

                    res += calculate_mad(result)
                    counter += 1
                    fitted_parameter.append([nOfCon, Dc, sAEC, result[0], result[1]])

        # Quasi-Random Walk
        for nOfCon in nOfCons:
            for Dc in [6000]:
                for sAEC in [1, 10]:
                    infection_scores, input_variables = retrieve_infection_scores_and_input_variables(input_file, cond_Dc=Dc, cond_sAEC=sAEC, cond_nOfCon=nOfCon)

                    if optimize:
                        result = minimize(objective_function_to_minimize, seed).x
                        print("qRW: nOfCon =", nOfCon, ", Dc =", Dc, ", sAEC =", sAEC, ":: Fitted parameter beta, gamma =", result[0], ",", result[1])
                    else:
                        result = [df[df["nOfCon"] == nOfCon][df["Dc"] == Dc][df["sAEC"] == sAEC]["beta"].values[0], df[df["nOfCon"] == nOfCon][df["Dc"] == Dc][df["sAEC"] == sAEC]["gamma"].values[0]]

                    res += calculate_mad(result)
                    counter += 1
                    fitted_parameter.append([nOfCon, Dc, sAEC, result[0], result[1]])

        # Save parameters to file
        if save_to_file:
            df = pd.DataFrame(fitted_parameter, columns=["nOfCon", "Dc", "sAEC", "beta", "gamma"])
            df.to_csv("fitted_parameter/CEF_"+system+".csv")

        print("CEF", system, "MAD =", res/counter)
