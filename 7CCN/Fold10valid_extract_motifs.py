# -*- coding: utf-8 -*-
"""
Created on Sun Jan  8 16:15:48 2023

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
output_file_path = "../outputs/valid10x10/motifs/"
directory = os.path.dirname(output_file_path)
if not os.path.exists(directory):
    print("Creating directory %s" % output_file_path)
    os.makedirs(directory)
else:
     print("Directory %s exists" % output_file_path)

batch_size = 100
num_filters = 300

# Device configuration
device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')

x = np.load('../data/one_hot_seqs.npy')
x = x.astype(np.float32)

y = np.load('../data/cell_type_array.npy')
y = y.astype(np.float32)

num_classes = y.shape[1]
peak_names = np.load('../data/peak_names.npy')

for model_index in (np.arange(10)+1):
    print('The running model is %s'%model_index)
    # load trained model
    model = aitac.ConvNet(num_classes, num_filters).to(device)
    checkpoint = torch.load('../outputs/valid10x10/' + 'model_%s'%model_index + '.ckpt')
    model.load_state_dict(checkpoint)
    
    #copy trained model weights to motif extraction model
    motif_model = aitac.motifCNN(model).to(device)
    motif_model.load_state_dict(model.state_dict())
    
    # Data loader
    dataset = torch.utils.data.TensorDataset(torch.from_numpy(x), torch.from_numpy(y))
    data_loader = torch.utils.data.DataLoader(dataset=dataset, batch_size=batch_size, shuffle=False)

    # run predictions with full model on all data
    pred_full_model, max_activations, activation_idx = aitac.test_model(data_loader, model, device, num_classes)
    correlations = plot_utils.plot_cors(y, pred_full_model, output_file_path, '_full_model_%s'%model_index)
    
    # find well predicted OCRs
    idx = np.argwhere(np.asarray(correlations)>0.75).squeeze()
    
    #get data subset for well predicted OCRs to run further test
    x2 = x[idx, :, :]
    y2 = y[idx, :]
    
    dataset = torch.utils.data.TensorDataset(torch.from_numpy(x2), torch.from_numpy(y2))
    data_loader = torch.utils.data.DataLoader(dataset=dataset, batch_size=batch_size, shuffle=False)
    
    # non-modified results for well-predicted OCRs only
    pred_full_model2 = pred_full_model[idx,:]
    correlations2 = plot_utils.plot_cors(y2, pred_full_model2, output_file_path, '_full_model_2_%s'%model_index)
    
    # get first layer activations and predictions with leave-one-filter-out
    activations, predictions = aitac.get_motifs(data_loader, motif_model, device, num_classes)
    
    np.save(output_file_path + "activations_model_%s.npy"%model_index, activations)
    np.save(output_file_path + "predictions_model_%s.npy"%model_index, predictions)
    np.save(output_file_path + "x2_model_%s.npy"%model_index, x2)
    np.save(output_file_path + "y2_model_%s.npy"%model_index, y2)
    np.save(output_file_path + "correlations2_model_%s.npy"%model_index, correlations2)
    np.save(output_file_path + "pred_full_model2_model_%s.npy"%model_index, pred_full_model2)
    
    # filt_corr, filt_infl, ave_filt_infl = plot_utils.plot_filt_corr(predictions, y2, correlations2, output_file_path)
    # infl, infl_by_OCR = plot_utils.plot_filt_infl(pred_full_model2, predictions, output_file_path)
    pwm, act_ind, nseqs, activated_OCRs, n_activated_OCRs, OCR_matrix = plot_utils.get_memes(activations, x2, y2, output_file_path, 'model_%s'%model_index)
    
    np.save(output_file_path + "pwm_model_%s.npy"%model_index, pwm)
    np.save(output_file_path + "act_ind_model_%s.npy"%model_index, act_ind)
    np.save(output_file_path + "nseqs_model_%s.npy"%model_index, nseqs)
    np.save(output_file_path + "activated_OCRs_model_%s.npy"%model_index, activated_OCRs)
    np.save(output_file_path + "n_activated_OCRs_model_%s.npy"%model_index, n_activated_OCRs)
    np.save(output_file_path + "OCR_matrix_model_%s.npy"%model_index, OCR_matrix)
