# -*- coding: utf-8 -*-
"""
Created on Sun Jan  8 16:02:24 2023

@author: Sonric
"""

import numpy as np
from numpy import random
import pandas as pd

from sklearn.model_selection import KFold
import torch.utils.data
import matplotlib
import os
import sys
matplotlib.use('Agg')

import aitac
import plot_utils

#create output directory
output_file_path = "../outputs/valid10x10/"
directory = os.path.dirname(output_file_path)
if not os.path.exists(directory):
    print("Creating directory %s" % output_file_path)
    os.makedirs(directory)
else:
     print("Directory %s exists" % output_file_path)

# Device configuration
device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')

# Hyper parameters
num_epochs = 20
batch_size = 100
learning_rate = 0.001
num_filters = 300

x = np.load('../data/one_hot_seqs.npy')
x = x.astype(np.float32)
y = np.load('../data/cell_type_array.npy')
y = y.astype(np.float32)
num_classes = y.shape[1]
peak_names = np.load('../data/peak_names.npy')


def cross_validate(x, y, peak_names, output_file_path):
    kf = KFold(n_splits=10, shuffle=True)

    pred_all = []
    corr_all = []
    peak_order = []
    
    model_index = 1
    for train_index, test_index in kf.split(x):
        
        print('The running model is %s'%model_index)
        
        train_data, eval_data = x[train_index, :, :], x[test_index, :, :]
        train_labels, eval_labels = y[train_index, :], y[test_index, :]
        train_names, eval_name = peak_names[train_index], peak_names[test_index]

        # Data loader
        train_dataset = torch.utils.data.TensorDataset(torch.from_numpy(train_data), torch.from_numpy(train_labels))
        train_loader = torch.utils.data.DataLoader(dataset=train_dataset, batch_size=batch_size, shuffle=False)

        eval_dataset = torch.utils.data.TensorDataset(torch.from_numpy(eval_data), torch.from_numpy(eval_labels))
        eval_loader = torch.utils.data.DataLoader(dataset=eval_dataset, batch_size=batch_size, shuffle=False)

        # create model 
        model = aitac.ConvNet(num_classes, num_filters).to(device)

        # Loss and optimizer
        criterion = aitac.pearson_loss
        optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)

        # train model
        model, best_loss = aitac.train_model(train_loader, eval_loader, model, device, criterion,  optimizer, num_epochs, output_file_path)

        # Predict on test set
        predictions, max_activations, max_act_index = aitac.test_model(eval_loader, model, device, num_classes)

        # plot the correlations histogram
        correlations = plot_utils.plot_cors(eval_labels, predictions, output_file_path, '_model_%s_'%model_index)

        pred_all.append(predictions)
        corr_all.append(correlations)
        peak_order.append(eval_name)
        
        torch.save(model.state_dict(), output_file_path + 'model_%s'%model_index + '.ckpt')
        torch.save(model, output_file_path + 'model_%s'%model_index + '.pth')
        model_index = model_index + 1
    
    pred_all = np.vstack(pred_all)
    corr_all = np.hstack(corr_all)
    peak_order = np.hstack(peak_order)

    return pred_all, corr_all, peak_order


predictions, correlations, peak_order = cross_validate(x, y, peak_names, output_file_path)

np.save(output_file_path + "predictions_trial.npy", predictions)
np.save(output_file_path + "correlations_trial.npy", correlations)
np.save(output_file_path + "peak_order.npy", peak_order)