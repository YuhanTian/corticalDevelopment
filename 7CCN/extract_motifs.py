import numpy as np
import torch
import torch.utils.data
import matplotlib
import os
import sys
matplotlib.use('Agg')

import aitac
import plot_utils

import time
from sklearn.model_selection import train_test_split

# Device configuration
device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')

# Hyper parameters
# num_classes = 23
batch_size = 100
num_filters = 300

#create output figure directory
model_name = 'test'
output_file_path = "../outputs/" + model_name + "/motifs/"
directory = os.path.dirname(output_file_path)
if not os.path.exists(directory):
    os.makedirs(directory)


# Load all data
x = np.load('../data/one_hot_seqs.npy')
x = x.astype(np.float32)

y = np.load('../data/cell_type_array.npy')
y = y.astype(np.float32)

num_classes = y.shape[1]
peak_names = np.load('../data/peak_names.npy')

# Data loader
dataset = torch.utils.data.TensorDataset(torch.from_numpy(x), torch.from_numpy(y))
data_loader = torch.utils.data.DataLoader(dataset=dataset, batch_size=batch_size, shuffle=False)

# load trained model
model = aitac.ConvNet(num_classes, num_filters).to(device)
checkpoint = torch.load('../models/' + model_name + '.ckpt')
model.load_state_dict(checkpoint)

#copy trained model weights to motif extraction model
motif_model = aitac.motifCNN(model).to(device)
motif_model.load_state_dict(model.state_dict())

# run predictions with full model on all data
pred_full_model, max_activations, activation_idx = aitac.test_model(data_loader, model, device, num_classes)
correlations = plot_utils.plot_cors(y, pred_full_model, output_file_path, 'full_model_')

# find well predicted OCRs
idx = np.argwhere(np.asarray(correlations)>0.75).squeeze()

#get data subset for well predicted OCRs to run further test
x2 = x[idx, :, :]
y2 = y[idx, :]

dataset = torch.utils.data.TensorDataset(torch.from_numpy(x2), torch.from_numpy(y2))
data_loader = torch.utils.data.DataLoader(dataset=dataset, batch_size=batch_size, shuffle=False)

# non-modified results for well-predicted OCRs only
pred_full_model2 = pred_full_model[idx,:]
correlations2 = plot_utils.plot_cors(y2, pred_full_model2, output_file_path, 'full_model_2')


# get first layer activations and predictions with leave-one-filter-out
start = time.time()
activations, predictions = aitac.get_motifs(data_loader, motif_model, device, num_classes)
print(time.time()- start)

np.save(output_file_path + "activations.npy", activations)
np.save(output_file_path + "predictions.npy", predictions)
np.save(output_file_path + "x2.npy", x2)
np.save(output_file_path + "y2.npy", y2)
np.save(output_file_path + "correlations2.npy", correlations2)
np.save(output_file_path + "pred_full_model2.npy", pred_full_model2)

activations = np.load(output_file_path + "activations.npy")
predictions = np.load(output_file_path + "predictions.npy")
x2 = np.load(output_file_path + "x2.npy")
y2 = np.load(output_file_path + "y2.npy")
correlations2 = np.load(output_file_path + "correlations2.npy")
pred_full_model2 = np.load(output_file_path + "pred_full_model2.npy")

filt_corr, filt_infl, ave_filt_infl = plot_utils.plot_filt_corr(predictions, y2, correlations2, output_file_path)

infl, infl_by_OCR = plot_utils.plot_filt_infl(pred_full_model2, predictions, output_file_path)

pwm, act_ind, nseqs, activated_OCRs, n_activated_OCRs, OCR_matrix = plot_utils.get_memes(activations, x2, y2, output_file_path)

np.save(output_file_path + "filt_corr.npy", filt_corr)
np.save(output_file_path + "filt_infl.npy", filt_infl)
np.save(output_file_path + "ave_filt_infl.npy", ave_filt_infl)
np.save(output_file_path + "infl.npy", infl)
np.save(output_file_path + "infl_by_OCR.npy", infl_by_OCR)
np.save(output_file_path + "pwm.npy", pwm)
np.save(output_file_path + "act_ind.npy", act_ind)
np.save(output_file_path + "nseqs.npy", nseqs)
np.save(output_file_path + "activated_OCRs.npy", activated_OCRs)
np.save(output_file_path + "n_activated_OCRs.npy", n_activated_OCRs)
np.save(output_file_path + "OCR_matrix.npy", OCR_matrix)

act_ind = np.load(output_file_path + "act_ind.npy")
activated_OCRs = np.load(output_file_path + "activated_OCRs.npy")
n_activated_OCRs = np.load(output_file_path + "n_activated_OCRs.npy")
OCR_matrix = np.load(output_file_path + "OCR_matrix.npy")
activations = np.load(output_file_path + "activations.npy")