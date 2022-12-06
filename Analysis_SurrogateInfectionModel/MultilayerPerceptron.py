#   Copyright by Christoph Saffer
#   Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
#   https://www.leibniz-hki.de/en/applied-systems-biology.html
#   HKI-Center for Systems Biology of Infection
#   Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
#   Adolf-Reichwein-Straße 23, 07745 Jena, Germany
#
#   This code is licensed under BSD 2-Clause
#   See the LICENSE file provided with this code for the full license.

import torch
import torch.nn as nn


# Class for feed forward Multilayer Perceptron (MLP) (Neural network)
class Feedforward(torch.nn.Module):
    def __init__(self, input_size, hidden_size, hidden_size2):
        super(Feedforward, self).__init__()
        self.input_size = input_size
        self.hidden_size = hidden_size
        self.hidden_size2 = hidden_size2
        self.hl1 = nn.Linear(self.input_size, self.hidden_size)
        self.hl2 = nn.Linear(self.hidden_size, self.hidden_size2)
        self.sigmoid = nn.Sigmoid()
        self.output_layer = torch.nn.Linear(self.hidden_size2, 1)

    def forward(self, x):
        hidden = self.hl1(x)
        sig = self.sigmoid(hidden)
        hidden2 = self.hl2(sig)
        sig2 = self.sigmoid(hidden2)
        output = self.output_layer(sig2)
        return torch.flatten(output)

    def optimize_net(self, data, target, epoch=25000, verbose=False):
        criterion = torch.nn.MSELoss()
        optimizer = torch.optim.Adam(self.parameters(), lr=0.0005)

        self.train()
        for epoch in range(epoch):
            optimizer.zero_grad()
            pred = self(data)  # get predictions NN(w)(x_i) of the current model on the data
            loss = criterion(pred.squeeze(), target)  # get MSE loss sum_i^n |NN(w)(x_i) - y_i|^2/n
            if verbose:
                if epoch % 1000 == 0:
                    print('Epoch {}: train loss: {}'.format(epoch, torch.sqrt(loss).item()))

            loss.backward()  # compute gradients
            optimizer.step()  # do a gradient descend step

    def rmse(self, data, target):
        criterion = torch.nn.MSELoss()
        pred = self(data)
        return torch.sqrt(criterion(pred.squeeze(), target)).item()

    def absloss(self, data, target):
        criterion = torch.nn.L1Loss()
        pred = self(data)
        return criterion(pred.squeeze(), target).item()
