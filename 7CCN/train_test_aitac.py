import numpy as np
from sklearn.model_selection import train_test_split
import torch.utils.data
import matplotlib
import os
import pathlib
import aitac

matplotlib.use('Agg')

# Device configuration
device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')

# Hyper parameters
num_epochs = 20
batch_size = 100
learning_rate = 0.001
num_filters = 300

#create output figure directory
# model_name = sys.argv[1]
model_name = 'test'
output_file_path = "../outputs/" + model_name + "/training/"

if not os.path.exists(output_file_path):
    print("Creating directory %s" % output_file_path)
    pathlib.Path(output_file_path).mkdir(parents=True, exist_ok=True) 
else:
     print("Directory %s exists" % output_file_path)

# load all data
x = np.load('../data/one_hot_seqs.npy')
x = x.astype(np.float32)
y = np.load('../data/cell_type_array.npy')
y = y.astype(np.float32)
num_classes = y.shape[1]
peak_names = np.load('../data/peak_names.npy')

# split the data into training and test sets
train_data, eval_data, train_labels, eval_labels, train_names, eval_names = train_test_split(x, y, peak_names, test_size=0.1, random_state=40)

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
model, best_loss = aitac.train_model(train_loader, eval_loader, model, device, criterion, optimizer, num_epochs, output_file_path)

# save the model checkpoint
models_file_path = "../models/"
if not os.path.exists(models_file_path):
    print("Creating directory %s" % models_file_path)
    os.makedirs(models_file_path)
else:
     print("Directory %s exists" % models_file_path)

torch.save(model.state_dict(), '../models/' + model_name + '.ckpt')

#save the whole model
torch.save(model, '../models/' + model_name + '.pth')

# Predict on test set
torch.load('../models/' + model_name + '.pth')

predictions_train, max_activations, max_act_index = aitac.test_model(train_loader, model, device, num_classes)
predictions_eval, max_activations, max_act_index = aitac.test_model(eval_loader, model, device, num_classes)

# #-------------------------------------------#
# #                 Save Files                #
# #-------------------------------------------#

np.save(output_file_path + "train_data.npy", train_data)
np.save(output_file_path + "eval_data.npy", eval_data)
np.save(output_file_path + "train_labels.npy", train_labels)
np.save(output_file_path + "eval_labels.npy", eval_labels)
np.save(output_file_path + "train_names.npy", train_names)
np.save(output_file_path + "eval_names.npy", eval_names)
np.save(output_file_path + "predictions_train.npy", predictions_train)
np.save(output_file_path + "predictions_eval.npy", predictions_eval)
