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
import glob


# Global variables
path_to_simulation_outcomes = ""  # path to configuration folder for hABM simulations, usually at configurations/namStudy/...
limit = 360  # 6 hours

# Change parameters to human or mouse
system = "human"  # or "mouse"
if system == "human":
    nOfMs = np.arange(2, 52, 2)
    nOfCons = [1, 2]
else:
    nOfMs = np.arange(0.1, 2.6, 0.1)
    nOfCons = [1, 2, 3]
sAECs = [1500, 5000, 15000, 50000, 150000, 500000]  # secretion rates
Dcs = [20, 60, 200, 600, 2000, 6000]  # diffusion coefficients


# Retrieves maximum simulation time for each run
def get_simulation_times(df: pd.DataFrame) -> list:
    simulation_durations = []
    maxRunId = max(set(df["runId"]))
    for i in range(1, maxRunId + 1):
        currMax = max(df[df["runId"] == i]["Time"])
        simulation_durations.append(currMax)
    return simulation_durations


# Calculates 95% confidence interval
def get_95_confidence_intverval(p: float, n: int) -> float:
    # https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval
    return 1.96 * np.sqrt((p * (1 - p)) / n)


# Calculates infection score and retrieves 95% confidence interval
def retrieve_IS_and_CI_from_simulation_output(file: str) -> list:
    df = pd.read_csv(file, names=["runId", "Time", "NumbCon", "NumbMacro"], sep=";", skiprows=[0], index_col=False)
    sim_times = get_simulation_times(df)
    inf_score = np.sum(np.array(sim_times) > limit) / len(sim_times)
    confidence_int = get_95_confidence_intverval(inf_score, len(sim_times))
    return [inf_score, confidence_int]


if __name__ == '__main__':

    if path_to_simulation_outcomes == "":
        print("ERROR: Path to simulation folder is not specified in global variables.")
    else:
        data = []
        # Loops over all configurations and reads measurement files from hABM
        # HIGHLY DEPENDS IN WHICH WAY THE DATA WAS SIMULATED
        for nOfCon in nOfCons:
            for nOfM in nOfMs:
                for Dc in Dcs:
                    for sAEC in sAECs:
                        try:
                            # ATTENTION: THE FOLLOWING LINE MUST BE ADJUSTED TO THE LOCAL FOLDERS
                            # IT READS THE FILE, BUT ONLY IF IT IS FOUND AT THE PROPOSED LOCATION
                            file = glob.glob(path_to_simulation_outcomes + "/*nOfM" + str(nOfM) + "_nOfCon" + str(nOfCon) + "_srAEC" + str(sAEC) + "_dc" + str(Dc) + "_*/measurements/*")[0]
                            df = pd.read_csv(file, names=["runId", "Time", "NumbCon", "NumbMacro"], sep=";", skiprows=[0], index_col=False)
                        except:
                            print("Can not read measurement file for nOfCon =", nOfCon, ", nOfM =", nOfM, ", Dc =", Dc, ", sAEC =", sAEC)
                            break

                        # Calculates infection scores and 95% confidence interval
                        infScore, infScoreCI = retrieve_IS_and_CI_from_simulation_output(file)
                        data.append([nOfCon, nOfM, Dc, sAEC, infScore, infScoreCI])

        # Generates pandas dataframe and writes it to csv file
        # simulated_data = pd.DataFrame(data, columns=["nOfCon", "nOfM", "Dc", "sAEC", "infScore", "confidence_int_95"])
        # simulated_data.to_csv("processed_data/"+system+"_processed_data.csv")
