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
import warnings
warnings.filterwarnings('ignore')  # setting ignore as a parameter


# Find the combinations of sAEC x Dc which have the same ratio sAEC/Dc
def find_indentical_ratios_sAEC_Dc(file: str) -> list:
    df = pd.read_csv(file, index_col=0)
    df["ratio"] = np.log(df["sAEC"] / df["Dc"])
    df = df.sort_values("ratio").reset_index()

    symm_values = []
    curr_values = []
    idx = 1
    for i in range(1, len(df)-1):
        last_ratio = df.loc[i - 1, "ratio"]
        curr_ratio = df.loc[i, "ratio"]
        if np.abs(last_ratio - curr_ratio) < 0.0001:
            if idx == 0:
                curr_values.append([df.loc[i - 1, "sAEC"], df.loc[i - 1, "Dc"]])
            curr_values.append([df.loc[i, "sAEC"], df.loc[i, "Dc"]])
            idx += 1
        else:
            idx = 0
            if len(curr_values) > 0:
                symm_values.append(curr_values)
            curr_values = []

    return symm_values


# Calculate the error between the infection scores of combinations sAEC x Dc with identical ratios sAEC/Dc
def calculate_symm_errors(values: list, processed_data: str) -> float:
    dists = []
    df = pd.read_csv(processed_data, index_col=0)
    nOfCons = df["nOfCon"].unique()
    nOfMs = df["nOfM"].unique()
    for nOfCon in nOfCons:
        ncons = []
        for nOfM in nOfMs:
            nofm = []
            for mult in values:
                for i in range(len(mult)):
                    for j in range(i + 1, len(mult)):
                        S1, D1 = mult[i]
                        S2, D2 = mult[j]
                        infS1 = df[df["nOfM"] == nOfM][df["nOfCon"] == nOfCon][df["sAEC"] == S1][df["Dc"] == D1]["infScore"].values[0]
                        infS2 = df[df["nOfM"] == nOfM][df["nOfCon"] == nOfCon][df["sAEC"] == S2][df["Dc"] == D2]["infScore"].values[0]
                        dist = np.abs(infS2 - infS1)
                        nofm.append(dist)
            ncons.append(np.mean(nofm))
        dists.append(ncons)

    return np.mean(dists)


if __name__ == "__main__":

    for system in ["human", "mouse"]:
        # Read in files
        input_file = "fitted_parameter/WSM_"+system+".csv"
        processed_data_file = "processed_data/"+system+"_processed_data.csv"

        # Find identical ratios
        idr = find_indentical_ratios_sAEC_Dc(input_file)

        # Print identical ratios
        # print("\nIdentical ratios", system, ": ")
        # for x in idr:
        #     print(x, " ~ log(sAEC/Dc) =", np.log(x[0][0]/x[0][1]))

        # Calculate symmetric error ME
        res = calculate_symm_errors(idr, processed_data_file)

        print("\nMean symmetric error", system, "ME =", res)

