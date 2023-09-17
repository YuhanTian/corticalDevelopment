# -*- coding: utf-8 -*-
"""
Created on Sat Dec 17 15:39:10 2022

@author: Sonric
"""

import pandas as pd
import preprocess_utils
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

multiomics_peak_region = pd.read_table('../data/peakRegion_summit251bp.bed', header=0)
multiomics_peak_region = multiomics_peak_region[['peakID','chr','start','end']]

multiomics_peak_height = pd.read_table('../data/peakRegion_cellTypeCounts.final.txt', header=0)

multiomics_peak_height = multiomics_peak_height[["chr","start","end","peakID","RG_Cyc","RG1","RG2","RG3","RG4","MigN","superficialCPN","Layer4CPN","deepCPN","CThPN","nIPC1","nIPC2","nIPC3","ImmatureN1","ImmatureN2","SCPN","SCPN_CSMN","SCPN_CTPN"]]
multiomics_cell_number = np.array([333,758,800,1092,1184,890,434,752,288,306,396,187,1048,1365,1929,2572,1044,797])

data_cell = np.log(multiomics_peak_height.iloc[:,4:]/multiomics_cell_number*np.sum(multiomics_cell_number)+1)

valid_peak_index = np.max(data_cell,axis=1) >= 4
data_cell_valid = data_cell.loc[valid_peak_index,:]
data_cell_valid_normalized,rank_mean = preprocess_utils.quantile_normalize(data_cell_valid.T)
data_cell_valid_normalized = data_cell_valid_normalized.T

np.save('../data/rank_mean.npy',rank_mean)

multiomics_peak_height_normalized = data_cell_valid_normalized
multiomics_peak_region_normalized = multiomics_peak_region.loc[valid_peak_index,:]


multiomics_peak_region_normalized.to_csv('../data/multiomics_peak_region_normalized.txt', sep='\t', index=False, header=False)

multiomics_peak_height_normalized = pd.concat([multiomics_peak_region_normalized['peakID'],multiomics_peak_height_normalized], axis=1)
multiomics_peak_height_normalized.to_csv('../data/multiomics_peak_height_normalized.csv', sep=',', index=False, header=False)