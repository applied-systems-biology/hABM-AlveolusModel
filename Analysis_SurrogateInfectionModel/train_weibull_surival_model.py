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


# Weibull survival model (WSM)
def weibull_survival_model(variables: list, parameters: float) -> float:
    nOfM, nOfCon, sAEC, Dc = variables
    lambda_p = parameters
    exponent = - lambda_p * nOfM

    return 1 - (1 - np.exp(exponent))**nOfCon


# Objective function to calibrate WSM parameters
def objective_function_to_minimize(parameters: np.array) -> float:
    list_of_squared_errors = ((infection_scores - weibull_survival_model(input_variables, parameters)) ** 2)
    return np.sum(list_of_squared_errors)/len(list_of_squared_errors)


# Calculate mean average deviation (MAD) for WSM
def calculate_mad(parameters: np.array) -> float:
    list_of_absolute_errors = np.abs(infection_scores - weibull_survival_model(input_variables, parameters))
    return np.sum(list_of_absolute_errors)/len(list_of_absolute_errors)


if __name__ == '__main__':

    for system in ["human", "mouse"]:
        # Set parameters for the corresponding system
        if system == "human":
            nOfMs = np.arange(2, 52, 2)  # human number of AM
            # nOfCons = [1, 2]  # including high fungal burden: human number of conidia
        else:
            nOfMs = np.arange(0.1, 2.6, 0.1)  # mouse number of AM
            # nOfCons = [1, 2, 3]  # including high fungal burden: mouse number of conidia
        nOfCons = [1]  # only low fungal burden
        sAECs = [1500, 5000, 15000, 50000, 150000, 500000]  # secretion rates
        Dcs = [20, 60, 200, 600, 2000, 6000]  # diffusion coefficients

        # Read in processed data
        input_file = "processed_data/" + system + "_processed_data.csv"  # path to retrieved infection scores
        res, counter, seed = 0, 0, [1]
        fitted_parameter = []

        # Read optimal parameters
        df = pd.read_csv("fitted_parameter/WSM_"+system+".csv", index_col=0)

        # Loop over all Dc and sAEC and calibrate WSM
        for Dc in Dcs:
            for sAEC in sAECs:
                if len(nOfCons) > 1:
                    infection_scores, input_variables = retrieve_infection_scores_and_input_variables(input_file, cond_Dc=Dc, cond_sAEC=sAEC)
                else:
                    infection_scores, input_variables = retrieve_infection_scores_and_input_variables(input_file, cond_Dc=Dc, cond_sAEC=sAEC, cond_nOfCon=1)

                if optimize:
                    result = minimize(objective_function_to_minimize, seed).x
                    print("Dc =", Dc, ", sAEC =", sAEC, ":: Fitted parameter Lambda =", result[0])
                else:
                    result = df[df["Dc"] == Dc][df["sAEC"] == sAEC]["Lambda"].values

                res += calculate_mad(result)
                counter += 1
                fitted_parameter.append([Dc, sAEC, result[0]])

        # Quasi-Random Walk
        for Dc in [6000]:
            for sAEC in [1, 10]:
                if len(nOfCons) > 1:
                    infection_scores, input_variables = retrieve_infection_scores_and_input_variables(input_file, cond_Dc=Dc, cond_sAEC=sAEC)
                else:
                    infection_scores, input_variables = retrieve_infection_scores_and_input_variables(input_file, cond_Dc=Dc, cond_sAEC=sAEC, cond_nOfCon=1)

                if optimize:
                    result = minimize(objective_function_to_minimize, seed).x
                    print("qRW: Dc =", Dc, ", sAEC =", sAEC, ":: Fitted parameter Lambda =", result[0])
                else:
                    result = df[df["Dc"] == Dc][df["sAEC"] == sAEC]["Lambda"].values

                res += calculate_mad(result)
                counter += 1
                fitted_parameter.append([Dc, sAEC, result[0]])

        # Save parameters to file
        if len(nOfCons) == 1 and save_to_file:
            df = pd.DataFrame(fitted_parameter, columns=["Dc", "sAEC", "Lambda"])
            df.to_csv("fitted_parameter/WSM_"+system+".csv")

        print("WSM", system, "MAD =", res/counter)
