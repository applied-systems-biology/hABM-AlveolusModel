#   Copyright by Christoph Saffer
#   Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
#   https://www.leibniz-hki.de/en/applied-systems-biology.html
#   HKI-Center for Systems Biology of Infection
#   Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
#   Adolf-Reichwein-Straße 23, 07745 Jena, Germany
#
#   This code is licensed under BSD 2-Clause
#   See the LICENSE file provided with this code for the full license.

from sklearn.ensemble import RandomForestRegressor
from MultilayerPerceptron import *
from train_surrogate_infection_model import *
from utils import retrieve_infection_scores_and_input_variables, retrieve_seed_point_for_SIM


# Global variables
save_to_file = False  # save parameters to file


# Make input variables and infection scores suitable for Multilayer Perceptron (MLP)
def to_tensor(variab: list, infsc: list) -> list:
    a = variab[0]
    b = np.log(variab[2]/variab[3])
    d = variab[1]
    e = infsc
    data = torch.tensor([a, b, d], dtype=torch.float32).t()
    target = torch.tensor(e, dtype=torch.float32)
    return data, target


# Make input variables suitable for Random Decision Forest (RDF)
def to_Scikit_array(variab: list) -> np.array:
    a = variab[0]
    b = np.log(variab[2]/variab[3])
    d = variab[1]
    return np.array([a, b, d]).transpose()


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
                # print("Run:", j*number_of_trials+(i+1), "/", number_of_second_trails*number_of_trials, ": Found new optimum for LSE =", optimal_value)

    return optimal_parameters, optimal_value


# Split dataset into train and test data
def validate(k=6, n=36, len_nOfM=25):

    # Divides data of length n * len_nOfM * nOfCons in k subsets
    # n: sAEC x Dc Combinations + qRW
    kfold_idx = []
    for nOfCon in nOfCons:
        numbers = list(np.arange(n))
        random.shuffle(numbers)
        numbers_test = np.array_split(numbers, k)
        kfold_idx.append(numbers_test)

    # Function that returns the idx-th set combination of train and test data for input variables and infection scores
    def split_test_train(idx):
        trainSetY, trainSetX, testSetY, testSetX = [], [], [], []

        for nOfCon in nOfCons:
            for i in range(n):
                if i in kfold_idx[nOfCon - 1][idx]:
                    for j in range(len_nOfM):
                        current_i = (nOfCon - 1) * n * len_nOfM + i * len_nOfM + j
                        testSetY.append(infection_scores0[current_i])
                        testSetX.append(input_variables0.T[current_i])
                else:
                    for j in range(len_nOfM):
                        current_i = (nOfCon - 1) * n * len_nOfM + i * len_nOfM + j
                        trainSetY.append(infection_scores0[current_i])
                        trainSetX.append(input_variables0.T[current_i])

        return np.array(trainSetY), np.array(trainSetX).transpose(), np.array(testSetY), np.array(testSetX).transpose()

    return split_test_train


if __name__ == '__main__':
    # Parameters for 5 times 6-fold cross validation
    number_of_sets = 6
    repeats = 5

    # Do cross validation for both systems
    for system in ["human", "mouse"]:
        if system == "human":
            nOfMs = np.arange(2, 52, 2)  # human number of AM
            nOfCons = [1, 2]  # including high fungal burden: human number of conidia
        else:
            nOfMs = np.arange(0.1, 2.6, 0.1)  # mouse number of AM
            nOfCons = [1, 2, 3]  # including high fungal burden: mouse number of conidi

        trerrors = []
        errors = []
        # Read in processed data
        input_file = "processed_data/" + system + "_processed_data.csv"  # path to retrieved infection scores
        # Retrieve infection scores and input variables from processed data (full dataset)
        infection_scores0, input_variables0 = retrieve_infection_scores_and_input_variables(input_file)
        # Retrieve seedpoints derived from generalized functions for SIM
        parameters = np.array(retrieve_seed_point_for_SIM(generalized_function_f, system))

        # Repeat cross validation 5 times
        for times in range(repeats):
            print("####", system, "#### Run", times + 1, "/", repeats)

            # Function split(i) ensures to provide the same splitted data for same index for each model
            split = validate(k=number_of_sets)

            # Train SIM on 5/6 of data and predict for 1/6 (test) - do this 6 times
            print("Train SIM")
            for i in range(number_of_sets):
                infScores, variables, infScoresT, variablesT = split(i)
                infection_scores, input_variables = infScores, variables
                par, val = train_surrogate_infection_model(parameters, number_of_trials=100, noise_to_seed=0.5)
                trerrors.append([calculate_mad(par), i, times, "train", "Surr", "mad"])
                infection_scores, input_variables = infScoresT, variablesT
                trerrors.append([calculate_mad(par), i, times, "test", "Surr", "mad"])

            # Train MLP on 5/6 of data and predict for 1/6 (test) - do this 6 times
            print("Train MLP")
            for i in range(number_of_sets):
                infs, v, infsT, vT = split(i)
                data, target = to_tensor(v, infs)
                data_test, target_test = to_tensor(vT, infsT)
                model = Feedforward(3, 13, 3)
                model.optimize_net(data, target, epoch=25000)

                train = model.absloss(data, target)
                test = model.absloss(data_test, target_test)
                trerrors.append([train, i, times, "train", "DNN", "mad"])
                trerrors.append([test, i, times, "test", "DNN", "mad"])

            # Train RDF on 5/6 of data and predict for 1/6 (test) - do this 6 times
            print("Train RDF")
            for i in range(number_of_sets):
                target, v, target_test, vT = split(i)
                data = to_Scikit_array(v)
                data_test = to_Scikit_array(vT)

                rf = RandomForestRegressor(n_estimators=100)
                rf.fit(data, target)

                train = np.sum(np.abs(rf.predict(data) - target)) / len(target)
                test = np.sum(np.abs(rf.predict(data_test) - target_test)) / len(target_test)
                trerrors.append([train, i, times, "train", "RF", "mad"])
                trerrors.append([test, i, times, "test", "RF", "mad"])

        # Save data to file
        if save_to_file:
            crossv = pd.DataFrame(trerrors, columns=["loss", "set", "run", "traintest", "model", "loss_type"])
            crossv.to_csv("cross_validation/cv_"+system+".csv")
